import { useState } from 'react'
import UploadZone from './components/UploadZone'
import DrugInput from './components/DrugInput'
import LoadingSpinner from './components/LoadingSpinner'
import ResultCard from './components/ResultCard'
import { analyzeVCF } from './services/api'

export default function App() {
  const [vcfFile, setVcfFile] = useState(null)
  const [drug, setDrug] = useState('')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState(null)
  const [error, setError] = useState(null)

  const handleAnalyze = async () => {
    if (!vcfFile) {
      setError('Please upload a VCF file first.')
      return
    }
    if (!drug.trim()) {
      setError('Please select or enter a drug name.')
      return
    }

    setError(null)
    setResult(null)
    setLoading(true)

    try {
      const data = await analyzeVCF(vcfFile, drug)
      setResult(data)
    } catch (err) {
      setError('Analysis failed. Please check the file and try again.')
    } finally {
      setLoading(false)
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

          {error && (
            <div className="error-box">{error}</div>
          )}

          <button
            className="analyze-btn"
            onClick={handleAnalyze}
            disabled={loading}
          >
            {loading ? 'Analyzing...' : 'Analyze Genetic Risk'}
          </button>
        </div>

        {loading && <LoadingSpinner />}
        {result && <ResultCard data={result} />}

      </div>
    </div>
  )
}