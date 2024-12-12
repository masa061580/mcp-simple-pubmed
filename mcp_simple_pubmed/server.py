"""
MCP server implementation for PubMed integration.
"""
import os
import json
import logging
import asyncio
from typing import Optional, Sequence, Dict, Any
from urllib.parse import urlparse, parse_qs

from mcp.server import Server
import mcp.types as types
from mcp.server.stdio import stdio_server
from .pubmed_search import PubMedSearch
from .pubmed_fetch import PubMedFetch

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-server")

app = Server("pubmed-server")

# Set up error handler
app.onerror = lambda error: logger.error(f"Server error: {error}")

# Initialize the clients
email = os.environ.get("PUBMED_EMAIL")
if not email:
    raise ValueError("PUBMED_EMAIL environment variable is required")
    
tool = os.environ.get("PUBMED_TOOL", "mcp-simple-pubmed")
api_key = os.environ.get("PUBMED_API_KEY")

pubmed_search = PubMedSearch(email=email, tool=tool, api_key=api_key)
pubmed_fetch = PubMedFetch()

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
- "breast cancer"[Title] AND "2024"[Date - Publication]
- "Smith J"[Author] AND "diabetes"
- "RNA"[MeSH Terms] AND "therapy"

Returns for each article:
- Paper titles
- Authors
- Publication details
- Abstract preview (when available)
- Web URLs for direct access
- DOI and PMC links when available
- URIs for accessing full text through this tool

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
            description="""Get full text and detailed information about a PubMed article using its ID.

This tool will attempt to:
1. Fetch the full text content (if available via PubMed Central)
2. Provide article metadata:
   - Title and abstract
   - Authors
   - Journal information
   - Publication date
   - Citation information
3. Generate access URLs:
   - PubMed web and mobile links
   - DOI link (if available)
   - PubMed Central link (if available)

If full text isn't directly available, the tool will provide information about 
where the paper can be accessed.

Example usage:
get_paper_fulltext(pmid="39661433")""",
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
            raise ValueError(f"Arguments must be a dictionary, got {type(arguments)}")

        if name == "search_pubmed":
            if "query" not in arguments:
                raise ValueError("Missing required argument: query")

            # Extract arguments
            query = arguments["query"]
            max_results = min(int(arguments.get("max_results", 10)), 50)

            # Perform the search
            results = await pubmed_search.search_articles(
                query=query,
                max_results=max_results
            )

            return [types.TextContent(
                type="text",
                text=json.dumps(results, indent=2)
            )]
            
        elif name == "get_paper_fulltext":
            if "pmid" not in arguments:
                raise ValueError("Missing required argument: pmid")
                
            # Get full text and metadata
            paper_info, urls = await pubmed_fetch.get_full_text(arguments["pmid"])
            
            # Combine results
            result = {
                "paper_info": paper_info,
                "urls": urls
            }
            
            return [types.TextContent(
                type="text",
                text=json.dumps(result, indent=2)
            )]
        
        else:
            raise ValueError(f"Unknown tool: {name}")
        
    except Exception as e:
        logger.exception(f"Error in call_tool: {str(e)}")
        return [types.TextContent(
            type="text",
            text=f"Error: {str(e)}",
            isError=True
        )]

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