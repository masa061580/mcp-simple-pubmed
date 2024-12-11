"""
Client for interacting with PubMed/Entrez API.
"""
import os
import time
import logging
import http.client
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Any, Tuple
from Bio import Entrez
from metapub import PubMedFetcher

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-client")

class PubMedClient:
    """Client for interacting with PubMed/Entrez API."""

    def __init__(self, email: str, tool: str, api_key: Optional[str] = None):
        """Initialize PubMed client with required credentials.

        Args:
            email: Valid email address for API access
            tool: Unique identifier for the tool
            api_key: Optional API key for higher rate limits
        """
        self.email = email
        self.tool = tool
        self.api_key = api_key
        
        # Configure Entrez
        Entrez.email = email
        Entrez.tool = tool
        if api_key:
            Entrez.api_key = api_key
            
        # Initialize metapub fetcher
        self.fetcher = PubMedFetcher()

    def _clean_text(self, text: Optional[str]) -> Optional[str]:
        """Clean and format text content."""
        if text is None:
            return None
        # Replace multiple spaces and newlines with single space
        cleaned = ' '.join(text.split())
        return cleaned

    async def search_articles(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """Search for articles matching the query.

        Args:
            query: Search query string
            max_results: Maximum number of results to return

        Returns:
            List of article metadata dictionaries
        """
        try:
            logger.info(f"Searching PubMed with query: {query}")
            results = []
            
            # Step 1: Search for article IDs
            handle = Entrez.esearch(db="pubmed", term=query, retmax=str(max_results))
            if not handle:
                logger.error("Got None handle from esearch")
                return []
                
            if isinstance(handle, http.client.HTTPResponse):
                logger.info("Got valid HTTP response from esearch")
                xml_content = handle.read()
                handle.close()
                
                # Parse XML to get IDs
                root = ET.fromstring(xml_content)
                id_list = root.findall('.//Id')
                
                if not id_list:
                    logger.info("No results found")
                    return []
                    
                pmids = [id_elem.text for id_elem in id_list]
                logger.info(f"Found {len(pmids)} articles")
                
                # Step 2: Get details for each article
                for pmid in pmids:
                    try:
                        article = self.fetcher.article_by_pmid(pmid)
                        results.append({
                            "pmid": pmid,
                            "title": article.title,
                            "abstract": article.abstract,
                            "journal": article.journal,
                            "authors": [str(author) for author in article.authors],
                            "publication_date": {
                                "year": article.year,
                                "month": article.month,
                                "day": article.day
                            },
                            "doi": article.doi,
                            "pmc_id": article.pmc,
                            "urls": self._generate_urls(pmid, article.doi, article.pmc),
                            "abstract_uri": f"pubmed://{pmid}/abstract",
                            "full_text_uri": f"pubmed://{pmid}/full_text"
                        })
                    except Exception as e:
                        logger.error(f"Error processing article {pmid}: {str(e)}")
                        continue
            
            return results

        except Exception as e:
            logger.exception(f"Error in search_articles: {str(e)}")
            raise
            
    async def get_full_text(self, pmid: str) -> Tuple[Dict[str, Any], Dict[str, str]]:
        """Get full text and metadata of an article using metapub.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            Tuple of (article info dict, URLs dict)
        """
        try:
            logger.info(f"Fetching article {pmid} using metapub")
            article = self.fetcher.article_by_pmid(pmid)
            urls = self._generate_urls(pmid, article.doi, article.pmc)
            
            result = {
                "pmid": pmid,
                "title": article.title,
                "journal": article.journal,
                "authors": [str(author) for author in article.authors],
                "year": article.year,
                "doi": article.doi,
                "pmc_id": article.pmc,
                "abstract": article.abstract
            }
            
            try:
                content = article.content
                if content:
                    result["full_text"] = self._clean_text(content)
                    result["content_source"] = "Full text available"
                else:
                    result["content_source"] = "Full text not directly available"
            except Exception as e:
                logger.warning(f"Could not get full text content: {str(e)}")
                result["content_source"] = f"Error getting full text: {str(e)}"
            
            # Add citation info if available
            if hasattr(article, 'citation'):
                result["citation"] = article.citation
                
            return result, urls
            
        except Exception as e:
            logger.exception(f"Error fetching article {pmid}: {str(e)}")
            return {
                "error": f"Error fetching article: {str(e)}",
                "pmid": pmid
            }, self._generate_urls(pmid)
            
    def _get_xml_text(self, elem: Optional[ET.Element], xpath: str) -> Optional[str]:
        """Helper method to safely get text from XML element."""
        if elem is None:
            return None
        found = elem.find(xpath)
        return found.text if found is not None else None
        
    def _generate_urls(self, pmid: str, doi: Optional[str] = None, pmc_id: Optional[str] = None) -> Dict[str, str]:
        """Generate URLs for human access.
        
        Args:
            pmid: PubMed ID
            doi: Optional DOI
            pmc_id: Optional PMC ID
            
        Returns:
            Dictionary with URLs
        """
        urls = {
            "pubmed": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            "pubmed_mobile": f"https://m.pubmed.ncbi.nlm.nih.gov/{pmid}/"
        }
        
        if doi:
            urls["doi"] = f"https://doi.org/{doi}"
        if pmc_id:
            urls["pmc"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"
            
        return urls