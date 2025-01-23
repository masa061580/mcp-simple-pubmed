"""
MCP server implementation for PubMed integration.
"""
import os
import json
import logging
from typing import Optional, Sequence, Dict, Any, Tuple
from urllib.parse import urlparse, parse_qs

from mcp.server import Server
import mcp.types as types
from mcp.server.stdio import stdio_server
from .pubmed_client import PubMedClient
from .fulltext_client import FullTextClient

import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stdin.reconfigure(encoding='utf-8')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-server")

app = Server("pubmed-server")

# Set up error handler
app.onerror = lambda error: logger.error(f"Server error: {error}")

def configure_clients() -> Tuple[PubMedClient, FullTextClient]:
    """Configure PubMed and full text clients with environment settings."""
    email = os.environ.get("PUBMED_EMAIL")
    if not email:
        raise ValueError("PUBMED_EMAIL environment variable is required")
        
    tool = os.environ.get("PUBMED_TOOL", "mcp-simple-pubmed")
    api_key = os.environ.get("PUBMED_API_KEY")

    pubmed_client = PubMedClient(email=email, tool=tool, api_key=api_key)
    fulltext_client = FullTextClient(email=email, tool=tool, api_key=api_key)
    
    return pubmed_client, fulltext_client

# Initialize the clients
pubmed_client, fulltext_client = configure_clients()

@app.list_tools()
async def list_tools() -> list[types.Tool]:
    """List available tools for interacting with PubMed."""
    return [
        types.Tool(
            name="search_pubmed",
            description="""Search PubMed for medical and life sciences research articles.

You can use these search features:
- Simple keyword search: "covid vaccine"
- Field-specific search:
  - Title search: [Title]
  - Author search: [Author]
  - MeSH terms: [MeSH Terms]
  - Journal: [Journal]
- Date ranges: Add year or date range like "2020:2024[Date - Publication]"
- Combine terms with AND, OR, NOT
- Use quotation marks for exact phrases

Examples:
- "covid vaccine" - basic search
- "breast cancer"[Title] AND "2023"[Date - Publication]
- "Smith J"[Author] AND "diabetes"
- "RNA"[MeSH Terms] AND "therapy"

The search will return:
- Paper titles
- Authors
- Publication details
- Abstract preview (when available)
- Links to full text (when available)
- DOI when available
- Keywords and MeSH terms

Note: Use quotes around multi-word terms for best results.""",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query to match against papers (e.g., 'covid vaccine', 'cancer treatment')"
                    },
                    "max_results": {
                        "type": "number",
                        "description": "Maximum number of results to return (default: 10)",
                        "default": 10,
                        "minimum": 1,
                        "maximum": 50
                    }
                },
                "required": ["query"]
            }
        ),
        types.Tool(
            name="get_paper_fulltext",
            description="""Get full text of a PubMed article using its ID.

        This tool attempts to retrieve the complete text of the paper if available through PubMed Central.
        If the paper is not available in PMC, it will return a message explaining why and provide information
        about where the text might be available (e.g., through DOI).

        Example usage:
        get_paper_fulltext(pmid="39661433")

        Returns:
        - If successful: The complete text of the paper
        - If not available: A clear message explaining why (e.g., "not in PMC", "requires journal access")""",
            inputSchema={
                "type": "object",
                "properties": {
                    "pmid": {
                        "type": "string",
                        "description": "PubMed ID of the article"
                    }
                },
                "required": ["pmid"]
            }
        )        
    ]

@app.call_tool()
async def call_tool(name: str, arguments: Dict[str, Any]) -> list[types.TextContent]:
    """Handle tool calls for PubMed operations."""
    try:
        # Log the received arguments for debugging
        logger.info(f"Received tool call: {name} with arguments: {json.dumps(arguments)}")
        
        # Validate arguments
        if not isinstance(arguments, dict):
            logger.error(f"Arguments must be a dictionary, got {type(arguments)}")
            return [types.TextContent(
                type="text",
                text="Invalid arguments: must be a dictionary",
                isError=True
            )]

        if name == "search_pubmed":
            # Validate search arguments
            if "query" not in arguments:
                logger.error("Missing required argument: query")
                return [types.TextContent(
                    type="text",
                    text="Missing required argument: query",
                    isError=True
                )]

            # Extract search arguments
            query = arguments["query"]
            max_results = min(int(arguments.get("max_results", 10)), 50)

            # Log the processed arguments
            logger.info(f"Processing search with query: {query}, max_results: {max_results}")

            # Perform the search
            results = await pubmed_client.search_articles(
                query=query,
                max_results=max_results
            )
            
            # Create resource URIs for articles
            articles_with_resources = []
            for article in results:
                pmid = article["pmid"]
                # Add original URIs
                article["abstract_uri"] = f"pubmed://{pmid}/abstract"
                article["full_text_uri"] = f"pubmed://{pmid}/full_text"
                
                # Add DOI URL if DOI exists
                if "doi" in article:
                    article["doi_url"] = f"https://doi.org/{article['doi']}"
                    
                # Add PubMed URLs
                article["pubmed_url"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                article["pubmed_fulltext_url"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmid}/"
                
                articles_with_resources.append(article)

            # Format the response
            formatted_results = json.dumps(articles_with_resources, indent=2)
            logger.info(f"Search completed successfully, found {len(results)} results")

            return [types.TextContent(
                type="text",
                text=formatted_results
            )]
            
        elif name == "get_paper_fulltext":
            # Validate fulltext arguments
            pmid = arguments.get("pmid")
            if not pmid:
                logger.error("Missing required argument: pmid")
                return [types.TextContent(
                    type="text",
                    text="Missing required argument: pmid",
                    isError=True
                )]

            # First check PMC availability
            available, pmc_id = await fulltext_client.check_full_text_availability(pmid)
            
            if available:
                full_text = await fulltext_client.get_full_text(pmid)
                if full_text:
                    logger.info(f"Successfully retrieved full text from PMC for PMID {pmid}")
                    return [types.TextContent(
                        type="text",
                        text=full_text
                    )]

            # Get article details to provide alternative locations
            article = await pubmed_client.get_article_details(pmid)
            
            message = "Full text is not available in PubMed Central.\n\n"
            message += "The article may be available at these locations:\n"
            message += f"- PubMed page: https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n"
            
            if article and "doi" in article:
                message += f"- Publisher's site (via DOI): https://doi.org/{article['doi']}\n"
                
            logger.info(f"Full text not available in PMC for PMID {pmid}, provided alternative locations")
            return [types.TextContent(
                type="text",
                text=message,
                isError=True
            )]
        
        else:
            logger.error(f"Unknown tool: {name}")
            return [types.TextContent(
                type="text",
                text=f"Unknown tool: {name}",
                isError=True
            )]
            
    except Exception as e:
        logger.exception(f"Error in call_tool ({name})")
        return [types.TextContent(
            type="text",
            text=f"Error processing request: {str(e)}",
            isError=True
        )]

@app.list_resources()
async def list_resources() -> list[types.Resource]:
    # For PubMed, resources are dynamic based on search results
    return []

@app.read_resource()
async def read_resource(uri: str) -> types.ReadResourceResult:
    try:
        parsed = urlparse(uri)
        if parsed.scheme != "pubmed":
            raise ValueError(f"Invalid URI scheme: {parsed.scheme}")

        # Extract PMID and resource type from URI
        parts = parsed.path.strip("/").split("/")
        if len(parts) != 2:
            raise ValueError(f"Invalid PubMed URI format: {uri}")

        pmid, resource_type = parts

        if resource_type == "abstract":
            article = await pubmed_client.get_article_details(pmid)
            return types.ReadResourceResult(
                contents=[types.ResourceContent(
                    uri=uri,
                    text=json.dumps(article, indent=2),
                    mimeType="application/json"
                )]
            )
        elif resource_type == "full_text":
            available, pmc_id = await fulltext_client.check_full_text_availability(pmid)
            if available:
                full_text = await fulltext_client.get_full_text(pmid)
                if full_text:
                    return types.ReadResourceResult(
                        contents=[types.ResourceContent(
                            uri=uri,
                            text=full_text,
                            mimeType="text/plain"
                        )]
                    )

            return types.ReadResourceResult(
                contents=[types.ResourceContent(
                    uri=uri,
                    text="Full text not available in PubMed Central",
                    mimeType="text/plain"
                )]
            )
        else:
            raise ValueError(f"Invalid resource type: {resource_type}")

    except Exception as e:
        logger.error(f"Error reading resource: {str(e)}")
        raise

async def main():
    """Run the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options()
        )

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    asyncio.run(main())