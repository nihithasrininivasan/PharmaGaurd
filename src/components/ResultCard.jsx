import { useState } from 'react'
import RiskBadge from './RiskBadge'

const SEVERITY_COLORS = {
  low:      { color: '#8C7B58', bg: 'rgba(140,123,88,0.15)',  label: 'LOW' },
  moderate: { color: '#DAA99A', bg: 'rgba(218,169,154,0.15)', label: 'MODERATE' },
  high:     { color: '#8D3437', bg: 'rgba(141,52,55,0.15)',   label: 'HIGH' },
  critical: { color: '#C75C5F', bg: 'rgba(199,92,95,0.2)',    label: 'CRITICAL' },
}

export default function ResultCard({ data }) {
  const [showExplain, setShowExplain] = useState(false)
  const severity = data.risk_assessment?.severity?.toLowerCase() || 'low'
  const sevStyle = SEVERITY_COLORS[severity] || SEVERITY_COLORS.low

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
          <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'flex-end', gap: 6 }}>
            <span style={{
              padding: '4px 12px',
              borderRadius: 100,
              fontSize: 10,
              fontFamily: 'var(--mono)',
              fontWeight: 700,
              letterSpacing: 2,
              color: sevStyle.color,
              border: `1.5px solid ${sevStyle.color}`,
              background: sevStyle.bg,
            }}>
              {sevStyle.label}
            </span>
            <div className="confidence-wrap">
              <div className="micro">Confidence</div>
              <div className="score">
                {Math.round((data.risk_assessment?.confidence_score || 0) * 100)}%
              </div>
            </div>
          </div>
        </div>

        <div className="gene-grid">
          <div className="micro">Pharmacogenomic Profile</div>
          <div className="gene-cols">
            <div className="gene-col">
              <div className="col-label">Gene</div>
              <div className="col-value">{data.pharmacogenomic_profile?.primary_gene}</div>
            </div>
            <div className="gene-col">
              <div className="col-label">Diplotype</div>
              <div className="col-value">{data.pharmacogenomic_profile?.diplotype}</div>
            </div>
            <div className="gene-col">
              <div className="col-label">Phenotype</div>
              <div className="col-value">{data.pharmacogenomic_profile?.phenotype}</div>
            </div>
          </div>
          {data.pharmacogenomic_profile?.detected_variants?.length > 0 && (
            <div style={{ marginTop: 14, borderTop: '1px solid var(--border)', paddingTop: 12 }}>
              <div className="col-label" style={{ marginBottom: 8 }}>Detected Variants</div>
              {data.pharmacogenomic_profile.detected_variants.map((v, i) => (
                <div key={i} style={{ fontSize: 12, color: 'var(--text-dim)', fontFamily: 'var(--mono)', marginBottom: 4 }}>
                  <span style={{ color: 'var(--accent-soft)' }}>{v.rsid}</span> — {v.effect}
                </div>
              ))}
            </div>
          )}
        </div>

        <div className="rec-box">
          <div className="rec-label">Clinical Recommendation</div>
          <p>{data.clinical_recommendation?.action}</p>
        </div>

        <button
          className="explain-toggle"
          onClick={() => setShowExplain(!showExplain)}
        >
          {showExplain ? '▲ Hide AI Explanation' : '▼ View AI Explanation'}
        </button>

        {showExplain && (
          <div className="explain-box">
            {data.llm_generated_explanation?.summary ||
              'AI-generated explanation will appear here once the model processes your genetic profile.'}
          </div>
        )}

        {data.quality_metrics && (
          <div style={{
            display: 'flex',
            alignItems: 'center',
            gap: 8,
            padding: '10px 14px',
            background: 'var(--bg2)',
            borderRadius: 8,
            border: '1px solid var(--border)',
          }}>
            <span style={{
              width: 8, height: 8, borderRadius: '50%',
              background: data.quality_metrics.vcf_parsing_success ? '#8C7B58' : '#C75C5F',
              flexShrink: 0,
            }} />
            <span style={{ fontSize: 11, fontFamily: 'var(--mono)', color: 'var(--text-dim)', letterSpacing: 1 }}>
              VCF PARSING: {data.quality_metrics.vcf_parsing_success ? 'PASS' : 'FAIL'}
            </span>
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