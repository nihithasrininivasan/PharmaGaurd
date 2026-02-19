import { useState, useRef, useEffect } from "react"
import { askPharmaGuard } from "../services/api"

const GLOW_COLORS = {
    critical: { shadow: "0 0 20px rgba(199,92,95,0.7), 0 0 40px rgba(199,92,95,0.3)", pulse: true, dot: "var(--red)" },
    high: { shadow: "0 0 20px rgba(199,92,95,0.5), 0 0 30px rgba(199,92,95,0.2)", pulse: true, dot: "var(--red)" },
    moderate: { shadow: "0 0 16px rgba(218,169,154,0.6), 0 0 30px rgba(218,169,154,0.2)", pulse: false, dot: "var(--accent)" },
    low: { shadow: "0 0 12px rgba(140,123,88,0.5), 0 0 24px rgba(140,123,88,0.2)", pulse: false, dot: "var(--green)" },
    none: { shadow: "0 0 12px rgba(140,123,88,0.3)", pulse: false, dot: "var(--green)" },
}

export default function AskPharmaGuard({ data }) {
    const [open, setOpen] = useState(false)
    const [question, setQuestion] = useState("")
    const [messages, setMessages] = useState([])
    const [loading, setLoading] = useState(false)

    const chatRef = useRef(null)

    useEffect(() => {
        if (chatRef.current) {
            chatRef.current.scrollTop = chatRef.current.scrollHeight
        }
    }, [messages, loading])

    if (!data) return null

    const severity = data.risk_assessment?.severity?.toLowerCase() || "none"
    const glow = GLOW_COLORS[severity] || GLOW_COLORS.none

    async function handleAsk() {
        if (!question.trim()) return

        const userMsg = { role: "user", text: question }
        setMessages(prev => [...prev, userMsg])
        setQuestion("")
        setLoading(true)

        const res = await askPharmaGuard({
            question: userMsg.text,
            gene: data.pharmacogenomic_profile?.primary_gene || "",
            diplotype: data.pharmacogenomic_profile?.diplotype || "",
            phenotype: data.pharmacogenomic_profile?.phenotype || "",
            drug: data.drug || "",
        })

        const aiMsg = { role: "ai", text: res.answer }
        setMessages(prev => [...prev, aiMsg])
        setLoading(false)
    }

    function handleKeyDown(e) {
        if (e.key === "Enter") handleAsk()
    }

    return (
        <div style={{
            position: "fixed",
            bottom: 24,
            right: 24,
            zIndex: 50,
            fontFamily: "var(--sans)",
        }}>

            {/* Floating Toggle Button */}
            <button
                onClick={() => setOpen(!open)}
                style={{
                    background: "var(--cyan)",
                    color: "var(--bg)",
                    padding: "12px 22px",
                    borderRadius: 100,
                    fontFamily: "var(--mono)",
                    fontSize: 12,
                    fontWeight: 700,
                    letterSpacing: 2,
                    textTransform: "uppercase",
                    border: "none",
                    cursor: "pointer",
                    boxShadow: glow.shadow,
                    animation: glow.pulse ? "pulse-glow 2s ease-in-out infinite" : "none",
                    transition: "all 0.3s ease",
                    display: "flex",
                    alignItems: "center",
                    gap: 8,
                    float: "right",
                }}
            >
                <span style={{ fontSize: 16 }}>💊</span>
                {open ? "Close" : "Ask PharmaGuard"}
            </button>

            {/* Chat Panel */}
            {open && (
                <div
                    className="animate-slideIn"
                    style={{
                        marginTop: 12,
                        width: 360,
                        background: "var(--card)",
                        border: "1px solid var(--border)",
                        borderRadius: 16,
                        overflow: "hidden",
                        boxShadow: "0 16px 48px rgba(0,0,0,0.5)",
                        float: "right",
                        clear: "both",
                        display: "flex",
                        flexDirection: "column",
                    }}
                >
                    {/* Header */}
                    <div style={{
                        background: "var(--card2)",
                        borderBottom: "1px solid var(--border)",
                        padding: "14px 18px",
                        display: "flex",
                        alignItems: "center",
                        gap: 10,
                        flexShrink: 0,
                    }}>
                        <span style={{
                            width: 8, height: 8, borderRadius: "50%",
                            background: glow.dot, flexShrink: 0,
                            animation: glow.pulse ? "pulse-glow 2s ease-in-out infinite" : "none",
                        }} />
                        <span style={{
                            fontFamily: "var(--mono)", fontSize: 11,
                            letterSpacing: 2, textTransform: "uppercase",
                            color: "var(--text-dim)",
                        }}>
                            Clinical AI Assistant
                        </span>
                    </div>

                    {/* Context Bar */}
                    <div style={{
                        padding: "10px 18px",
                        borderBottom: "1px solid var(--border)",
                        display: "flex", flexWrap: "wrap", gap: 6,
                        flexShrink: 0,
                    }}>
                        {data.pharmacogenomic_profile?.primary_gene && (
                            <span style={{
                                padding: "2px 8px", borderRadius: 100, fontSize: 10,
                                fontFamily: "var(--mono)", fontWeight: 700,
                                background: "rgba(176,124,212,0.15)", color: "#b07cd4",
                            }}>
                                {data.pharmacogenomic_profile.primary_gene}
                            </span>
                        )}
                        {data.drug && (
                            <span style={{
                                padding: "2px 8px", borderRadius: 100, fontSize: 10,
                                fontFamily: "var(--mono)", fontWeight: 700,
                                background: "rgba(140,123,88,0.15)", color: "#8C7B58",
                            }}>
                                {data.drug}
                            </span>
                        )}
                        {data.pharmacogenomic_profile?.phenotype && (
                            <span style={{
                                padding: "2px 8px", borderRadius: 100, fontSize: 10,
                                fontFamily: "var(--mono)", fontWeight: 700,
                                background: "rgba(218,169,154,0.15)", color: "#DAA99A",
                            }}>
                                {data.pharmacogenomic_profile.phenotype}
                            </span>
                        )}
                    </div>

                    {/* Scrollable Chat Area */}
                    <div
                        ref={chatRef}
                        style={{
                            height: 260,
                            overflowY: "auto",
                            padding: 14,
                            display: "flex",
                            flexDirection: "column",
                            gap: 10,
                        }}
                    >
                        {messages.length === 0 && !loading && (
                            <div style={{
                                padding: "24px 12px",
                                textAlign: "center",
                                color: "var(--text-dim)",
                                fontSize: 12,
                                fontFamily: "var(--mono)",
                                opacity: 0.6,
                            }}>
                                Ask any clinical question about the patient's pharmacogenomic profile.
                            </div>
                        )}

                        {messages.map((m, i) => (
                            <div
                                key={i}
                                style={{
                                    alignSelf: m.role === "user" ? "flex-end" : "flex-start",
                                    maxWidth: "85%",
                                    padding: "10px 14px",
                                    borderRadius: m.role === "user" ? "14px 14px 4px 14px" : "14px 14px 14px 4px",
                                    background: m.role === "user" ? "var(--cyan)" : "var(--bg2)",
                                    color: m.role === "user" ? "var(--bg)" : "var(--text-dim)",
                                    fontSize: 13,
                                    lineHeight: 1.6,
                                    border: m.role === "user" ? "none" : "1px solid var(--border)",
                                    animation: "slideIn 0.2s ease-out",
                                }}
                            >
                                {m.text}
                            </div>
                        ))}

                        {loading && (
                            <div style={{
                                alignSelf: "flex-start",
                                padding: "10px 14px",
                                borderRadius: "14px 14px 14px 4px",
                                background: "var(--bg2)",
                                border: "1px solid var(--border)",
                                fontSize: 12,
                                color: "var(--text-dim)",
                                fontFamily: "var(--mono)",
                                letterSpacing: 1,
                            }}>
                                Thinking clinically...
                            </div>
                        )}
                    </div>

                    {/* Input Area */}
                    <div style={{
                        borderTop: "1px solid var(--border)",
                        padding: 12,
                        display: "flex",
                        gap: 8,
                        flexShrink: 0,
                    }}>
                        <input
                            style={{
                                flex: 1,
                                background: "var(--bg2)",
                                border: "1.5px solid var(--border)",
                                borderRadius: 10,
                                padding: "10px 14px",
                                fontSize: 13,
                                color: "var(--text)",
                                fontFamily: "var(--sans)",
                                outline: "none",
                                transition: "border-color 0.2s",
                                boxSizing: "border-box",
                            }}
                            placeholder="Ask a clinical question..."
                            value={question}
                            onChange={(e) => setQuestion(e.target.value)}
                            onKeyDown={handleKeyDown}
                            onFocus={(e) => e.target.style.borderColor = "var(--cyan)"}
                            onBlur={(e) => e.target.style.borderColor = "var(--border)"}
                        />
                        <button
                            onClick={handleAsk}
                            disabled={loading || !question.trim()}
                            style={{
                                background: loading ? "var(--border2)" : "var(--cyan)",
                                color: loading ? "var(--text-dim)" : "var(--bg)",
                                fontFamily: "var(--mono)",
                                fontSize: 11,
                                fontWeight: 700,
                                letterSpacing: 1,
                                borderRadius: 10,
                                padding: "10px 16px",
                                border: "none",
                                cursor: loading ? "not-allowed" : "pointer",
                                transition: "all 0.2s",
                                flexShrink: 0,
                            }}
                        >
                            Ask
                        </button>
                    </div>

                    {/* Disclaimer */}
                    <div style={{
                        borderTop: "1px solid var(--border)",
                        padding: "8px 18px",
                        fontSize: 10,
                        color: "var(--text-dim)",
                        fontStyle: "italic",
                        opacity: 0.7,
                        lineHeight: 1.4,
                        flexShrink: 0,
                    }}>
                        AI-assisted insights only. Consult a licensed clinician for decisions.
                    </div>
                </div>
            )}
        </div>
    )
}
