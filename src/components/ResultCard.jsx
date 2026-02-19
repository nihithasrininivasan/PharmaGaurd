import { useState } from 'react'
import RiskBadge from './RiskBadge'

export default function ResultCard({ data }) {
  const [showExplain, setShowExplain] = useState(false)

  return (
    <div className="result-card">
      <div className="result-header">
        <div className="result-header-left">
          <div className="label">Patient ID</div>
          <div className="value">{data.patient_id}</div>
        </div>

        <div className="result-header-right">
          <div className="label">Drug</div>
          <div className="drug-name">{data.drug}</div>
        </div>
      </div>

      <div className="result-body">
        <div className="risk-row">
          <div className="risk-label-wrap">
            <div className="micro">Risk Assessment</div>
            <RiskBadge riskLabel={data.risk_assessment?.risk_label} />
          </div>

          <div className="confidence-wrap">
            <div className="micro">Confidence</div>
            <div className="score">
              {Math.round((data.risk_assessment?.confidence_score || 0) * 100)}%
            </div>
          </div>
        </div>

        <div className="gene-grid">
          <div className="micro">Pharmacogenomic Profile</div>
          <div className="gene-cols">
            <div className="gene-col">
              <div className="col-label">Gene</div>
              <div className="col-value">
                {data.pharmacogenomic_profile?.primary_gene}
              </div>
            </div>

            <div className="gene-col">
              <div className="col-label">Diplotype</div>
              <div className="col-value">
                {data.pharmacogenomic_profile?.diplotype}
              </div>
            </div>

            <div className="gene-col">
              <div className="col-label">Phenotype</div>
              <div className="col-value">
                {data.pharmacogenomic_profile?.phenotype}
              </div>
            </div>
          </div>
        </div>

        <div className="rec-box">
          <div className="rec-label">Clinical Recommendation</div>
          <p>{data.clinical_recommendation?.action}</p>
        </div>

        <button
          className="explain-toggle"
          onClick={() => setShowExplain(!showExplain)}
        >
          {showExplain ? 'Hide Explanation' : 'View AI Explanation'}
        </button>

        {showExplain && (
          <div className="explain-box">
            {data.llm_generated_explanation?.summary}
          </div>
        )}

        <div className="disclaimer">
          This tool provides AI-assisted pharmacogenomic insights.
          Final prescribing decisions must be made by licensed clinicians.
        </div>
      </div>
    </div>
  )
}
