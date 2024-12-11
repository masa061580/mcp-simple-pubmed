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

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-server")

app = Server("pubmed-server")

# Set up error handler
app.onerror = lambda error: logger.error(f"Server error: {error}")

def configure_pubmed_client() -> PubMedSearch:
    """Configure PubMed client with environment settings."""
    email = os.environ.get("PUBMED_EMAIL")
    if not email:
        raise ValueError("PUBMED_EMAIL environment variable is required")
        
    tool = os.environ.get("PUBMED_TOOL", "mcp-simple-pubmed")
    api_key = os.environ.get("PUBMED_API_KEY")

    return PubMedSearch(email=email, tool=tool, api_key=api_key)

# Initialize the client
pubmed_client = configure_pubmed_client()

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
            results = await pubmed_client.search_articles(
                query=query,
                max_results=max_results
            )

            return [types.TextContent(
                type="text",
                text=json.dumps(results, indent=2)
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