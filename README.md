# MCP Simple PubMed

An MCP server that provides access to PubMed articles through the Entrez API.

## Features

- Search PubMed database using keywords
- Access article abstracts
- Download full text when available (for open access articles) *not implemented yet!*

## Installation

```bash
pip install mcp-simple-pubmed
```

## Configuration

The server requires the following environment variables:

- `PUBMED_EMAIL`: Your email address (required by NCBI)
- `PUBMED_API_KEY`: Optional API key for higher rate limits 

The standard rate limit is 3 requests / second. No rate limiting was implemented, as it is highly unlikely in the typical usage scenario that your AI would generate more traffic. If you need it, you can [register for an API key](https://www.ncbi.nlm.nih.gov/account/) which will give you 10 requests / second. Read about [this on NCBI pages](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen).

## Usage with Claude Desktop

Add to your Claude Desktop configuration (`claude_desktop_config.json`):

(Mac OS)

```json
{
  "mcpServers": {
    "simple-pubmed": {
      "command": "python",
      "args": ["-m", "mcp_simple_pubmed"],
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "PUBMED_API_KEY": "your-api-key" 
      }
    }
  }
}
```

(Windows)


```json
{
  "mcpServers": {
    "simple-pubmed": {
      "command": "C:\\Users\\YOUR_USERNAME\\AppData\\Local\\Programs\\Python\\Python311\\python.exe",
      "args": [
        "-m",
        "mcp_simple_pubmed"
      ],
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "PUBMED_API_KEY": "your-api-key" 
      }
    }
  }
}
```

## License

MIT License