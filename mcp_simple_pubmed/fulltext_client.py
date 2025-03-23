"""
Client for retrieving full text content of PubMed articles.
Handles both XML retrieval and text extraction with robust error handling.
"""
import logging
import time
from typing import Optional, Tuple, Dict, Any
from Bio import Entrez, Medline
import xml.etree.ElementTree as ET

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-fulltext")

class FullTextClient:
    """Client for retrieving full text content from PubMed Central."""

    def __init__(self, email: str, tool: str, api_key: Optional[str] = None):
        """Initialize full text client with required credentials."""
        self.email = email
        self.tool = tool
        self.api_key = api_key
        
        # Configure Entrez
        Entrez.email = email
        Entrez.tool = tool
        if api_key:
            Entrez.api_key = api_key

    async def check_full_text_availability(self, pmid: str) -> Tuple[bool, Optional[str]]:
        """Check if full text is available in PMC and get PMC ID if it exists."""
        try:
            logger.info(f"Checking PMC availability for PMID {pmid}")
            
            # First try to get article details to check for direct PMC ID
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            try:
                record = Medline.read(handle)
                # Check if PMC ID is directly in the record
                if 'PMC' in record:
                    pmc_id = record['PMC']
                    logger.info(f"Found PMC ID {pmc_id} directly in record for PMID {pmid}")
                    return True, pmc_id
            except Exception as e:
                logger.warning(f"Error reading Medline record: {str(e)}")
            finally:
                handle.close()
            
            # If not found in record, try elink
            link_handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
            link_results = Entrez.read(link_handle)
            link_handle.close()
            
            # Parse result to find PMC ID
            if link_results and len(link_results) > 0:
                linksets = link_results[0]
                if 'LinkSetDb' in linksets and len(linksets['LinkSetDb']) > 0:
                    links = linksets['LinkSetDb'][0]
                    if 'Link' in links and len(links['Link']) > 0:
                        pmc_id = links['Link'][0]['Id']
                        logger.info(f"Found PMC ID {pmc_id} via elink for PMID {pmid}")
                        return True, pmc_id
            
            logger.info(f"No PMC ID found for PMID {pmid}")
            return False, None
            
        except Exception as e:
            logger.exception(f"Error checking PMC availability for PMID {pmid}: {str(e)}")
            return False, None

    async def get_full_text(self, pmid: str, format_as_text: bool = False) -> Optional[str]:
        """Get full text of the article if available through PMC."""
        try:
            # First check availability and get PMC ID
            available, pmc_id = await self.check_full_text_availability(pmid)
            if not available or pmc_id is None:
                logger.info(f"Full text not available in PMC for PMID {pmid}")
                return None

            logger.info(f"Fetching full text for PMC ID {pmc_id}")
            
            # Use a single request for handling documents
            full_text_handle = Entrez.efetch(
                db="pmc", 
                id=pmc_id, 
                rettype="xml",
                retmode="xml"
            )
            
            content = full_text_handle.read()
            full_text_handle.close()
            
            # Convert bytes to string if needed
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            
            # Format as text if requested
            if format_as_text:
                try:
                    return self._extract_text_from_pmc_xml(content)
                except Exception as text_error:
                    logger.error(f"Error extracting text from XML: {str(text_error)}")
                    # Return the XML as fallback
                    return content
            
            return content
            
        except Exception as e:
            logger.exception(f"Error getting full text for PMID {pmid}: {str(e)}")
            return None
    
    def _extract_text_from_pmc_xml(self, xml_content: str) -> str:
        """Extract readable text content from PMC XML."""
        try:
            root = ET.fromstring(xml_content)
            
            # Dictionary for text parts
            parts = {}
            
            # Get article title
            title_elem = root.find(".//article-title")
            if title_elem is not None and title_elem.text:
                parts['title'] = self._clean_text(title_elem.text)
            
            # Get abstract and body text...
            # [Implementation details...]
            
            return "\n\n".join(text_parts)
            
        except ET.ParseError as e:
            logger.error(f"Error parsing PMC XML: {str(e)}")
            raise ValueError(f"Could not parse PMC XML content: {str(e)}")