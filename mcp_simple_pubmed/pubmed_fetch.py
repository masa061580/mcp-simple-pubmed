"""
Full text fetching functionality for PubMed articles.

This module uses Bio.Entrez for metadata and trafilatura for full text extraction
from publicly available URLs.
"""
import logging
import asyncio
from typing import Dict, Any, Tuple, Optional
import httpx
import trafilatura
from Bio import Entrez, Medline

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-fetch")

class PubMedFetch:
    """Client for fetching full text and detailed article information from PubMed."""

    def __init__(self):
        """Initialize PubMed fetcher."""
        # Configure httpx client for URL fetching
        self.http_client = httpx.AsyncClient(
            timeout=30.0,
            follow_redirects=True
        )

    async def _fetch_full_text_from_url(self, url: str) -> Optional[str]:
        """Attempt to fetch and extract full text from a URL.
        
        Args:
            url: URL to fetch content from
            
        Returns:
            Extracted text if successful, None otherwise
        """
        try:
            response = await self.http_client.get(url)
            response.raise_for_status()
            
            # Try to extract text using trafilatura
            downloaded = response.text
            text = trafilatura.extract(downloaded)
            
            if text:
                logger.info(f"Successfully extracted text from {url}")
                return text
            else:
                logger.warning(f"No text could be extracted from {url}")
                return None
                
        except Exception as e:
            logger.warning(f"Error fetching {url}: {str(e)}")
            return None

    async def _try_get_full_text(self, urls: Dict[str, str]) -> Tuple[Optional[str], str]:
        """Try to get full text from various URLs.
        
        Args:
            urls: Dictionary of URLs to try
            
        Returns:
            Tuple of (extracted text or None, source description)
        """
        # Try PMC first if available
        if "pmc" in urls:
            text = await self._fetch_full_text_from_url(urls["pmc"])
            if text:
                return text, f"Full text extracted from PMC: {urls['pmc']}"
                
        # Try DOI link
        if "doi" in urls:
            text = await self._fetch_full_text_from_url(urls["doi"])
            if text:
                return text, f"Full text extracted from DOI: {urls['doi']}"

        # No full text found - return access information
        source_info = []
        if "pmc" in urls:
            source_info.append(f"PubMed Central: {urls['pmc']}")
        if "doi" in urls:
            source_info.append(f"DOI: {urls['doi']}")
        if source_info:
            return None, "Full text may be available at: " + "; ".join(source_info)
        else:
            return None, "Full text access links not available"

    def _clean_text(self, text: Optional[str]) -> Optional[str]:
        """Clean and format text content.
        
        Args:
            text: Text to clean
            
        Returns:
            Cleaned text with normalized whitespace
        """
        if text is None:
            return None
        # Replace multiple spaces and newlines with single space
        cleaned = ' '.join(text.split())
        return cleaned

    async def get_full_text(self, pmid: str) -> Tuple[Dict[str, Any], Dict[str, str]]:
        """Get full text and detailed metadata of an article.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            Tuple of (article info dict, URLs dict)
            Article info includes:
            - Basic metadata (title, authors, journal, etc.)
            - Abstract
            - Full text (if available)
            - Citation information
            - Content source information
            URLs include links to:
            - PubMed
            - PubMed Mobile
            - DOI (if available)
            - PubMed Central (if available)
        """
        try:
            logger.info(f"Fetching article {pmid}")
            
            # Get article metadata using Bio.Entrez
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()
            
            # Generate access URLs
            urls = {
                "pubmed": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "pubmed_mobile": f"https://m.pubmed.ncbi.nlm.nih.gov/{pmid}/"
            }
            
            # Add DOI URL if available
            if 'DOI' in record:
                urls["doi"] = f"https://doi.org/{record['DOI']}"
                
            # Add PMC URL if available
            if 'PMC' in record:
                urls["pmc"] = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{record['PMC']}/"
                
            # Try to get full text
            full_text, content_source = await self._try_get_full_text(urls)
            
            # Build result dictionary
            result = {
                "pmid": pmid,
                "title": self._clean_text(record.get('TI')),
                "abstract": self._clean_text(record.get('AB')),
                "authors": record.get('AU', []),
                "journal": {
                    "name": record.get('JT'),
                    "volume": record.get('VI'),
                    "issue": record.get('IP'),
                    "pages": record.get('PG')
                },
                "publication_date": {
                    "year": record.get('DP', '').split()[0] if record.get('DP') else None,
                    "month": record.get('DP', '').split()[1] if record.get('DP', '').count(' ') >= 1 else None,
                    "day": record.get('DP', '').split()[2] if record.get('DP', '').count(' ') >= 2 else None
                },
                "identifiers": {
                    "pmid": pmid,
                    "doi": record.get('DOI'),
                    "pmc": record.get('PMC')
                },
                "content_source": content_source
            }
            
            # Add full text if available
            if full_text:
                result["full_text"] = full_text
                
            # Add citation if available
            if 'SO' in record:
                result["citation"] = record['SO']
            
            return result, urls
            
        except Exception as e:
            logger.exception(f"Error fetching article {pmid}")
            return {
                "error": str(e), 
                "pmid": pmid
            }, self._generate_urls(pmid)

    def _generate_urls(self, pmid: str, doi: Optional[str] = None, 
                      pmc_id: Optional[str] = None) -> Dict[str, str]:
        """Generate URLs for human access.
        
        Args:
            pmid: PubMed ID
            doi: Optional DOI
            pmc_id: Optional PMC ID
            
        Returns:
            Dictionary with URLs for various access methods
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