import axios from 'axios'

const BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

const MOCK_RESPONSE = {
  patient_id: 'PATIENT_001',
  drug: 'CODEINE',
  timestamp: new Date().toISOString(),
  risk_assessment: {
    risk_label: 'Toxic',
    confidence_score: 0.91,
    severity: 'critical',
  },
  pharmacogenomic_profile: {
    primary_gene: 'CYP2D6',
    diplotype: '*1/*2xN',
    phenotype: 'URM',
    detected_variants: [
      { rsid: 'rs1080985', effect: 'Ultra-Rapid Metabolizer allele' }
    ],
  },
  clinical_recommendation: {
    action: 'Avoid Codeine. Switch to a non-CYP2D6-metabolized opioid such as morphine or oxymorphone at standard doses.',
    cpic_level: 'Strong',
  },
  llm_generated_explanation: {
    summary: 'This patient carries duplicated CYP2D6*2 alleles, making them an Ultra-Rapid Metabolizer. Codeine is converted to morphine by CYP2D6 and in this genotype the conversion is accelerated, producing toxic morphine blood levels. This causes a high risk of respiratory depression. CPIC strongly recommends avoiding codeine entirely in Ultra-Rapid Metabolizers.',
  },
  quality_metrics: {
    vcf_parsing_success: true,
  },
}

const USE_MOCK = true

export async function analyzeVCF(vcfFile, drugs) {
  if (USE_MOCK) {
    await new Promise((r) => setTimeout(r, 2200))
    return { ...MOCK_RESPONSE, drug: drugs }
  }

  const formData = new FormData()
  formData.append('file', vcfFile)
  formData.append('drugs', drugs)

  const response = await axios.post(`${BASE_URL}/analyze`, formData, {
    headers: { 'Content-Type': 'multipart/form-data' },
  })
  return response.data
}