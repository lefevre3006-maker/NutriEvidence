import streamlit as st
import os
from groq import Groq
from Bio import Entrez

# --- CONFIGURATION ---
st.set_page_config(page_title="Nutri-Evidence GPT", page_icon="🥗", layout="centered")

# Hide standard Streamlit style for better embedding
hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_st_style, unsafe_allow_html=True)

# --- SIDEBAR (API KEY) ---
# For a public tool on your site, you usually provide the key yourself via Secrets
# so the user doesn't have to enter one.
if "GROQ_API_KEY" in st.secrets:
    api_key = st.secrets["GROQ_API_KEY"]
    client = Groq(api_key=api_key)
else:
    st.error("System Error: API Key missing.")
    st.stop()

# Email for PubMed (Required)
Entrez.email = "nutrition.services@bread-and-better.com"

# --- LOGIC ---
def get_evidence(query):
    # 1. Search PubMed (Last 5 Years, Clinical Trials)
    # We add "Clinical Trial" to the query automatically
    search_term = f"{query} AND (Clinical Trial[ptyp] OR Randomized Controlled Trial[ptyp]) AND 2020:2025[dp]"
    
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list:
            return None, "No recent clinical trials found on PubMed for this specific combination."

        # 2. Fetch Abstracts
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles = Entrez.read(handle)
        
        context_text = ""
        sources = []
        
        for article in articles['PubmedArticle']:
            try:
                title = article['MedlineCitation']['Article']['ArticleTitle']
                # Get abstract safely
                abstract_list = article['MedlineCitation']['Article']['Abstract']['AbstractText']
                abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
                
                context_text += f"\n- STUDY: {title}\n  ABSTRACT: {abstract}\n"
                
                # Store source for display
                pmid = article['MedlineCitation']['PMID']
                sources.append(f"[{title}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
            except:
                continue
                
        return sources, context_text

    except Exception as e:
        return None, f"Search Error: {str(e)}"

# --- UI LAYOUT ---
st.image("https://cdn-icons-png.flaticon.com/512/3004/3004458.png", width=80)
st.title("Nutri-Evidence GPT")
st.markdown("**The Scientific Fact-Checker.** Enter a supplement or food to analyze the latest Clinical Trials (2020-2025).")

user_query = st.text_input("What do you want to check?", placeholder="e.g., Curcumin for Arthritis, Vitamin D for Sleep...")

if st.button("Analyze Science"):
    if not user_query:
        st.warning("Please enter a query.")
    else:
        with st.spinner("Scanning PubMed & Reading Abstracts..."):
            sources, context = get_evidence(user_query)
            
            if sources:
                # AI Analysis
                system_prompt = """
                You are a Senior Scientist. Synthesize the provided clinical trial abstracts.
                1. What is the consensus? (Positive/Inconclusive/Negative)
                2. Highlight strength of evidence (Sample sizes, p-values).
                3. Be objective. Do not sell anything.
                """
                
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
                
                with st.expander("📚 View Sources (PubMed)"):
                    for s in sources:
                        st.markdown(f"- {s}")
            else:
                st.error(context)

st.divider()
st.caption("⚠️ **Disclaimer:** This tool uses Artificial Intelligence to summarize public scientific data. It is for informational purposes only and does not constitute medical advice. Always consult a doctor.")