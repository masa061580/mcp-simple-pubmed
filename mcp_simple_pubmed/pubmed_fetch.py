"""
Full text fetching functionality for PubMed articles.

This module focuses solely on retrieving full text content from PMC
using Bio.Entrez. It does not handle metadata retrieval which is
done by pubmed_search.py.
"""
import logging
from typing import Optional
import xml.etree.ElementTree as ET
from Bio import Entrez, Medline

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed-fetch")

class PubMedFetch:
    """Client for fetching full text from PubMed Central."""

    def _extract_text_from_pmc_xml(self, xml_content: str) -> str:
        """Extract readable text content from PMC XML.
        
        Args:
            xml_content: PMC article XML
            
        Returns:
            Extracted text content
        """
        try:
            root = ET.fromstring(xml_content)
            
            # Get article title
            title = root.find(".//article-title")
            text_parts = []
            if title is not None:
                text_parts.append(title.text)
            
            # Get abstract
            abstract = root.find(".//abstract")
            if abstract is not None:
                for p in abstract.findall(".//p"):
                    if p.text:
                        text_parts.append(p.text)
            
            # Get main body text
            body = root.find(".//body")
            if body is not None:
                for p in body.findall(".//p"):
                    if p.text:
                        text_parts.append(p.text)
                        
            return "\n\n".join(text_parts)
            
        except ET.ParseError as e:
            logger.error(f"Error parsing PMC XML: {str(e)}")
            raise ValueError(f"Could not parse PMC XML content: {str(e)}")
        except Exception as e:
            logger.error(f"Error extracting text from PMC XML: {str(e)}")
            raise ValueError(f"Error processing PMC content: {str(e)}")

    async def get_full_text(self, pmid: str) -> str:
        """Get full text of an article if available.
        
        Args:
            pmid: PubMed ID of the article
            
        Returns:
            Full text content if available, otherwise an error message
            explaining why the text is not available.
            
        Raises:
            ValueError: If there are issues accessing or parsing the content
        """
        try:
            # First get PMC ID if available
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
            record = Medline.read(handle)
            handle.close()
            
            if 'PMC' in record:
                logger.info(f"Found PMC ID for article {pmid}, attempting to fetch full text")
                pmc_id = record['PMC']
                
                # Get full text from PMC
                pmc_handle = Entrez.efetch(db='pmc', id=pmc_id, rettype='full', retmode='xml')
                full_text = pmc_handle.read()
                pmc_handle.close()
                
                # Parse XML and extract text
                text = self._extract_text_from_pmc_xml(full_text)
                if text:
                    return text
                else:
                    return "Error: No text content found in PMC XML"
                    
            elif 'DOI' in record:
                return f"Full text not available in PMC. Article has DOI {record['DOI']} - full text may be available through publisher"
            else:
                return "Full text not available - article is not in PMC and has no DOI"
                
        except Exception as e:
            logger.exception(f"Error getting full text for article {pmid}")
            return f"Error retrieving full text: {str(e)}"