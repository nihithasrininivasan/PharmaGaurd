from pydantic import BaseModel, Field
from typing import List, Optional

class VariantCall(BaseModel):
    chrom: str
    pos: int
    rsid: Optional[str] = None
    ref: str
    alt: str
    zygosity: str = Field(..., description="HET, HOM_REF, or HOM_ALT")
    quality: float
    filter: Optional[str] = None

class GenotypeData(BaseModel):
    sample_id: str
    gene_symbol: str
    variants: List[VariantCall]
    coverage_mean: Optional[float] = None
