import { useState } from 'react'
import UploadZone from './components/UploadZone'
import DrugInput from './components/DrugInput'
import LoadingSpinner from './components/LoadingSpinner'
import ResultCard from './components/ResultCard'
import { analyzeVCF } from './services/api'

const LOADING_STEPS = [
  'Uploading VCF...',
  'Analyzing Pharmacogenomics...',
  'Generating Clinical Explanation...',
]

export default function App() {
  const [vcfFile, setVcfFile] = useState(null)
  const [drug, setDrug] = useState('')
  const [loading, setLoading] = useState(false)
  const [loadingStep, setLoadingStep] = useState(0)
  const [result, setResult] = useState(null)
  const [jobId, setJobId] = useState(null)
  const [error, setError] = useState(null)

  const handleAnalyze = async () => {
    if (!vcfFile) { setError('Please upload a VCF file first.'); return }
    if (!drug.trim()) { setError('Please select or enter a drug name.'); return }

    setError(null)
    setResult(null)
    setLoading(true)
    setLoadingStep(0)

    const stepTimer1 = setTimeout(() => setLoadingStep(1), 800)
    const stepTimer2 = setTimeout(() => setLoadingStep(2), 1800)

    try {
      const { data, jobId: jid } = await analyzeVCF(vcfFile, drug)
      setResult(data)
      setJobId(jid)
    } catch (err) {
      setError('Analysis failed. Please check the file and try again.')
    } finally {
      clearTimeout(stepTimer1)
      clearTimeout(stepTimer2)
      setLoading(false)
      setLoadingStep(0)
    }
  }

  return (
    <div className="app">
      <div className="container">

        <div className="header">
          <p className="header-logo">RIFT 2026 / Pharmacogenomics Track</p>
          <h1>Pharma<span>Guard</span></h1>
          <p>Genetic risk prediction for safer prescribing</p>
        </div>

        <div className="card">
          <UploadZone onFileSelected={setVcfFile} />
          <DrugInput value={drug} onChange={setDrug} />

          {error && <div className="error-box">{error}</div>}

          <button
            className="analyze-btn"
            onClick={handleAnalyze}
            disabled={loading}
          >
            {loading ? LOADING_STEPS[loadingStep] : 'Analyze Genetic Risk'}
          </button>
        </div>

        {loading && <LoadingSpinner message={LOADING_STEPS[loadingStep]} />}
        {result && <ResultCard data={result} jobId={jobId} />}

      </div>
    </div>
  )
}