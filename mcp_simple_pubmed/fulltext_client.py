"""
Client for retrieving full text content of PubMed articles.
Separate from main PubMed client to maintain code separation and stability.
"""
import logging
import xml.etree.ElementTree as ET
import http.client
from typing import Optional
from Bio import Entrez

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-fulltext")

class FullTextClient:
    """Client for retrieving full text content from PubMed Central."""

    def __init__(self, email: str, tool: str, api_key: Optional[str] = None):
        """Initialize full text client with required credentials.

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

    async def get_pmcid(self, pmid: str) -> Optional[str]:
        """Get PMC ID for a given PubMed ID if available.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            PMC ID if available, None otherwise
        """
        try:
            logger.info(f"Looking up PMC ID for PMID {pmid}")
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
            if not handle:
                logger.info(f"No PMC link found for PMID {pmid}")
                return None
                
            xml_content = handle.read()
            handle.close()
            
            # Parse XML to get PMC ID
            root = ET.fromstring(xml_content)
            pmc_id_elem = root.find('.//LinkSetDb/Link/Id')
            if pmc_id_elem is None:
                logger.info(f"No PMC ID found for PMID {pmid}")
                return None
                
            pmc_id = pmc_id_elem.text
            logger.info(f"Found PMC ID {pmc_id} for PMID {pmid}")
            return pmc_id
            
        except Exception as e:
            logger.exception(f"Error getting PMC ID for PMID {pmid}: {str(e)}")
            return None

    async def get_full_text(self, pmid: str) -> Optional[str]:
        """Get full text of the article if available through PMC.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            Full text content if available, None otherwise
        """
        try:
            # First get the PMC ID
            pmc_id = await self.get_pmcid(pmid)
            if not pmc_id:
                logger.info(f"No PMC ID available for PMID {pmid}")
                return None
            
            logger.info(f"Fetching full text for PMC ID {pmc_id}")
            full_text_handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="text")
            if not full_text_handle:
                logger.info(f"Could not fetch full text for PMC ID {pmc_id}")
                return None
                
            full_text = full_text_handle.read()
            full_text_handle.close()
            
            if isinstance(full_text, bytes):
                full_text = full_text.decode('utf-8')
                
            return full_text
            
        except Exception as e:
            logger.exception(f"Error getting full text for PMID {pmid}: {str(e)}")
            return None

    async def check_full_text_availability(self, pmid: str) -> dict:
        """Check if full text is available and get access information.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            Dictionary with availability information:
            {
                "has_pmc": boolean,
                "pmc_id": string or None,
                "status": "available" | "not_in_pmc" | "error",
                "message": string
            }
        """
        try:
            pmc_id = await self.get_pmcid(pmid)
            if not pmc_id:
                return {
                    "has_pmc": False,
                    "pmc_id": None,
                    "status": "not_in_pmc",
                    "message": "Article not available in PubMed Central"
                }
                
            return {
                "has_pmc": True,
                "pmc_id": pmc_id,
                "status": "available",
                "message": "Full text available in PubMed Central"
            }
            
        except Exception as e:
            logger.exception(f"Error checking full text availability for PMID {pmid}: {str(e)}")
            return {
                "has_pmc": False,
                "pmc_id": None,
                "status": "error",
                "message": f"Error checking availability: {str(e)}"
            }
