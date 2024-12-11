# MCP Simple PubMed

An MCP server that provides access to PubMed articles through the Entrez API.

## Features

- Search PubMed database using keywords
- Access article abstracts
- Download full text when available (for open access articles)
- Rate limiting compliance with NCBI guidelines

## Installation

```bash
pip install mcp-simple-pubmed
```

## Configuration

The server requires the following environment variables:

- `PUBMED_EMAIL`: Your email address (required by NCBI)
- `PUBMED_TOOL`: Tool identifier (defaults to 'mcp-simple-pubmed')
- `PUBMED_API_KEY`: Optional API key for higher rate limits

## Usage with Claude Desktop

Add to your Claude Desktop configuration (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "pubmed": {
      "command": "mcp-simple-pubmed",
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "PUBMED_API_KEY": "your-api-key"  // optional
      }
    }
  }
}
```

## Development

To develop or modify this server:

1. Clone the repository
2. Install dependencies: `pip install -e .`
3. Run the server: `python -m mcp_simple_pubmed`

## License

MIT License