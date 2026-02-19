import axios from 'axios'

const BASE_URL = 'http://localhost:8000'


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
    summary: 'This patient carries duplicated <gene>CYP2D6</gene> <diplotype>*2xN</diplotype> alleles, making them an Ultra-Rapid Metabolizer. <drug>Codeine</drug> is converted to morphine by <gene>CYP2D6</gene> and in this genotype the conversion is accelerated, producing toxic morphine blood levels. The variant <variant>rs1080985</variant> is responsible for this effect. <cpic>CPIC</cpic> strongly recommends avoiding <drug>Codeine</drug> entirely in Ultra-Rapid Metabolizers.',
  },
  quality_metrics: {
    vcf_parsing_success: true,
    extra_metadata: {
      heatmap_intensity: 3
    }
  },
}

const USE_MOCK = false

export async function analyzeVCF(vcfFile, drugs) {
  if (USE_MOCK) {
    await new Promise((r) => setTimeout(r, 2200))
    return { data: { ...MOCK_RESPONSE, drug: drugs }, jobId: null }
  }

  const formData = new FormData()
  formData.append('drug', drugs)
  formData.append('vcf', vcfFile)
  formData.append('patient_id', 'PATIENT_001')

  const response = await axios.post(`${BASE_URL}/api/v1/analyze`, formData, {
    headers: { 'Content-Type': 'multipart/form-data' },
  })

  const data = response.data

  // Extract job_id from placeholder summary if present
  const jobMatch = data.llm_generated_explanation?.summary?.match(/job_id:(.*)/)
  const jobId = jobMatch ? jobMatch[1].trim() : null

  return { data, jobId }
}

export async function fetchExplanation(jobId) {
  const response = await axios.get(`${BASE_URL}/api/v1/explanation/${jobId}`)
  return response.data.summary
}