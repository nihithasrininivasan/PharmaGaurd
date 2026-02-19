import { useCallback, useState } from 'react'
import { useDropzone } from 'react-dropzone'

export default function UploadZone({ onFileSelected }) {
  const [fileName, setFileName] = useState(null)
  const [error, setError] = useState(null)

  const onDrop = useCallback((acceptedFiles, rejectedFiles) => {
    setError(null)
    if (rejectedFiles.length > 0) {
      setError('Invalid file. Please upload a .vcf file under 5MB.')
      return
    }
    const file = acceptedFiles[0]
    setFileName(file.name)
    onFileSelected(file)
  }, [onFileSelected])

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: { 'text/plain': ['.vcf'] },
    maxSize: 5 * 1024 * 1024,
    multiple: false,
  })

  return (
    <div>
      <p className="section-label">Genetic Data File</p>
      <div
        {...getRootProps()}
        className={
          'upload-zone' +
          (isDragActive ? ' active' : '') +
          (fileName ? ' has-file' : '')
        }
      >
        <input {...getInputProps()} />
        <div className="upload-icon">
          <svg viewBox="0 0 24 24">
            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
            <polyline points="17 8 12 3 7 8" />
            <line x1="12" y1="3" x2="12" y2="15" />
          </svg>
        </div>
        {fileName ? (
          <>
            <p className="file-name">{fileName}</p>
            <p className="sub">File loaded. Ready to analyze.</p>
          </>
        ) : (
          <>
            <p>{isDragActive ? 'Drop your VCF file here' : 'Drag and drop your .vcf file here'}</p>
            <p className="sub">or click to browse - max 5MB</p>
          </>
        )}
      </div>
      {error && <p className="upload-error">{error}</p>}
    </div>
  )
}