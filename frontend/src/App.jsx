import { useState } from 'react'
import UploadZone from './components/UploadZone'
import DrugInput from './components/DrugInput'
import LoadingSpinner from './components/LoadingSpinner'
import ResultCard from './components/ResultCard'
import AskPharmaGuard from './components/AskPharmaGuard'
import { analyzeVCF } from './services/api'

const LOADING_STEPS = [
  'Uploading VCF...',
  'Analyzing Pharmacogenomics...',
  'Generating Clinical Explanation...',
]

export default function App() {
  const [vcfFile, setVcfFile] = useState(null)
  const [selectedDrugs, setSelectedDrugs] = useState([])
  const [loading, setLoading] = useState(false)
  const [loadingStep, setLoadingStep] = useState(0)
  const [results, setResults] = useState([])
  const [error, setError] = useState(null)

  const toggleDrug = (drug) => {
    setSelectedDrugs((prev) =>
      prev.includes(drug)
        ? prev.filter((d) => d !== drug)
        : [...prev, drug]
    )
  }

  const handleAnalyze = async () => {
    if (!vcfFile) { setError('Please upload a VCF file first.'); return }
    if (selectedDrugs.length === 0) { setError('Please select at least one drug.'); return }

    setError(null)
    setResults([])
    setLoading(true)
    setLoadingStep(0)

    const stepTimer1 = setTimeout(() => setLoadingStep(1), 800)
    const stepTimer2 = setTimeout(() => setLoadingStep(2), 1800)

    try {
      // Analyze each drug in parallel
      const promises = selectedDrugs.map(async (drug) => {
        try {
          const { data, jobId } = await analyzeVCF(vcfFile, drug)
          return { data, jobId, drug, error: null }
        } catch (err) {
          return {
            data: null,
            jobId: null,
            drug,
            error: err.response?.data?.detail || `Analysis failed for ${drug}.`,
          }
        }
      })

      const allResults = await Promise.all(promises)
      const successes = allResults.filter((r) => r.data)
      const failures = allResults.filter((r) => r.error)

      setResults(successes)

      if (failures.length > 0 && successes.length === 0) {
        setError(failures.map((f) => f.error).join('\n'))
      } else if (failures.length > 0) {
        setError(`Some analyses failed: ${failures.map((f) => f.drug).join(', ')}`)
      }
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
          <DrugInput selectedDrugs={selectedDrugs} onToggle={toggleDrug} />

          {error && <div className="error-box">{error}</div>}

          <button
            className="analyze-btn"
            onClick={handleAnalyze}
            disabled={loading}
          >
            {loading
              ? LOADING_STEPS[loadingStep]
              : `Analyze Genetic Risk${selectedDrugs.length > 1 ? ` (${selectedDrugs.length} drugs)` : ''}`}
          </button>
        </div>

        {loading && <LoadingSpinner message={LOADING_STEPS[loadingStep]} />}
        {results.map((r, i) => (
          <ResultCard key={r.drug + i} data={r.data} jobId={r.jobId} />
        ))}
        {results.length > 0 && <AskPharmaGuard data={results[0].data} />}

      </div>
    </div>
  )
}