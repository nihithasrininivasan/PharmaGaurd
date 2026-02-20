const SUPPORTED_DRUGS = [
  'CODEINE',
  'WARFARIN',
  'CLOPIDOGREL',
  'SIMVASTATIN',
  'AZATHIOPRINE',
  'FLUOROURACIL',
]

export default function DrugInput({ selectedDrugs, onToggle }) {
  return (
    <div className="drug-section">
      <p className="section-label">Drug Selection <span style={{ fontSize: 11, color: 'var(--text-dim)', fontWeight: 400 }}>(select one or more)</span></p>
      <div className="drug-chips">
        {SUPPORTED_DRUGS.map((drug) => (
          <button
            key={drug}
            onClick={() => onToggle(drug)}
            className={'drug-chip' + (selectedDrugs.includes(drug) ? ' selected' : '')}
          >
            {drug}
          </button>
        ))}
      </div>
    </div>
  )
}