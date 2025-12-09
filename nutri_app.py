import streamlit as st
import os
from groq import Groq
from Bio import Entrez

# --- CONFIGURATION ---
st.set_page_config(page_title="Nutri-Evidence GPT", page_icon="🥗", layout="centered")

# Hide standard Streamlit style for clean embedding
hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown("""
<style>
    .main {
        background-color: #B2C9AD;
    }
    .stApp {
        background-color: #fbffe4 !important;
    }
</style>
""", unsafe_allow_html=True)

# --- SIDEBAR (API KEY) ---
if "GROQ_API_KEY" in st.secrets:
    api_key = st.secrets["GROQ_API_KEY"]
    client = Groq(api_key=api_key)
else:
    st.error("System Error: API Key missing.")
    st.stop()

Entrez.email = "bioexpertise.contact@gmail.com"

# --- LOGIC ---
def get_evidence(substance, benefit):
    # 1. CONSTRUCT QUERY
    # We combine inputs and look for high-evidence paper types
    base_query = f"{substance}"
    if benefit:
        base_query += f" AND {benefit}"
    
    # Broader Search: 2015-2025, Reviews, Meta-Analyses, and Trials
    search_term = f"{base_query} AND (Clinical Trial[ptyp] OR Meta-Analysis[ptyp] OR Review[ptyp] OR Randomized Controlled Trial[ptyp]) AND 2015:2025[dp]"
    
    try:
        # 2. SEARCH PUBMED (Fetch top 10 most relevant)
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=10, sort="relevance")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list:
            return None, "No relevant scientific papers found in the last 10 years."

        # 3. FETCH ABSTRACTS
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles = Entrez.read(handle)
        
        context_text = ""
        sources = []
        
        for article in articles['PubmedArticle']:
            try:
                citation = article['MedlineCitation']
                title = citation['Article']['ArticleTitle']
                pub_year = citation['Article']['Journal']['JournalIssue']['PubDate'].get('Year', 'N/A')
                
                # Extract Abstract safely
                abs_list = citation['Article']['Abstract']['AbstractText']
                abstract = " ".join(abs_list) if isinstance(abs_list, list) else str(abs_list)
                
                context_text += f"\n- STUDY ({pub_year}): {title}\n  ABSTRACT: {abstract}\n"
                
                pmid = citation['PMID']
                sources.append(f"[{title} ({pub_year})](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
            except:
                continue
                
        return sources, context_text

    except Exception as e:
        return None, f"Search Error: {str(e)}"

# --- UI LAYOUT ---
st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=80)
st.title("Nutri-Evidence GPT")
st.markdown("**The Scientific Fact-Checker.** Analysis based on Clinical Trials & Meta-Analyses (2015-2025).")

# TWO COLUMN LAYOUT
col1, col2 = st.columns(2)

with col1:
    substance_input = st.text_input("1. Ingredient / Food", placeholder="e.g. Collagen")

with col2:
    benefit_input = st.text_input("2. Benefit / Condition (Optional)", placeholder="e.g. Hair Loss")

if st.button("Analyze Science"):
    if not substance_input:
        st.warning("Please enter at least an Ingredient.")
    else:
        with st.spinner(f"Analyzing scientific literature for '{substance_input}'..."):
            sources, context = get_evidence(substance_input, benefit_input)
            
            if sources:
                # AI Analysis
                system_prompt = """
                You are a Senior Scientist. Synthesize the provided abstracts from PubMed.
                
                Structure your answer:
                1. **Executive Summary:** Does it work? (Yes / No / Inconclusive).
                2. **Quality of Evidence:** Are these large meta-analyses or small trials? 
                3. **Detailed Findings:** Briefly explain the biological mechanism or specific results.
                
                Be objective. If the evidence is weak, say so clearly.
                """
                
                try:
                    completion = client.chat.completions.create(
                        messages=[
                            {"role": "system", "content": system_prompt},
                            {"role": "user", "content": context}
                        ],
                        model="llama-3.3-70b-versatile",
                        temperature=0.1
                    )
                    
                    st.success("Analysis Complete")
                    st.markdown("### 🧬 Scientific Verdict")
                    st.write(completion.choices[0].message.content)
                    
                    with st.expander(f"📚 Analyzed {len(sources)} Papers (Source Links)"):
                        for s in sources:
                            st.markdown(f"- {s}")
                            
                except Exception as e:
                     st.error(f"AI Analysis Error: {e}")

            else:
                st.error(context)

st.divider()
st.caption("⚠️ **Disclaimer:** This tool uses Artificial Intelligence to summarize public scientific data (PubMed). It is for informational purposes only and does not constitute medical advice.")
