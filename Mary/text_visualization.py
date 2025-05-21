{\rtf1\ansi\ansicpg1252\cocoartf2821
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 text_visualization.py\
\
# Using wordcloud and matplotlib\
from wordcloud import WordCloud\
import matplotlib.pyplot as plt\
from collections import Counter\
import numpy as np\
\
import seaborn as sns\
import networkx as nx\
from matplotlib_venn import venn2\
import pandas as pd\
\
# From scikit-learn\
from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS\
\
from gensim.models import LdaModel\
from gensim.corpora import Dictionary\
\
import numpy as np\
from sklearn.feature_extraction.text import TfidfVectorizer\
import seaborn as sns\
import re\
sklearn_stop_words = set(ENGLISH_STOP_WORDS)\
\
\
def compare_topics(text1, text2, num_topics=10, text1_name = "Topics in Text 1", text2_name = "Topics in Text 2"):\
    # Prepare texts\
    texts = [text1.lower().split(), text2.lower().split()]\
    dictionary = Dictionary(texts)\
    \
    # Create corpus\
    corpus = [dictionary.doc2bow(text) for text in texts]\
    \
    # Train model\
    lda = LdaModel(corpus, num_topics=num_topics, id2word=dictionary)\
    \
    # Get topics for each text\
    topics1 = lda[corpus[0]]\
    topics2 = lda[corpus[1]]\
    \
    # Plot comparison\
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))\
    \
    # Topic distribution for text1\
    t1_topics = sorted(topics1, key=lambda x: x[1], reverse=True)\
    ax[0].barh([dictionary[id] for id, _ in t1_topics], \
              [prob for _, prob in t1_topics])\
    ax[0].set_title(text1_name)\
    \
    # Topic distribution for text2\
    t2_topics = sorted(topics2, key=lambda x: x[1], reverse=True)\
    ax[1].barh([dictionary[id] for id, _ in t2_topics], \
              [prob for _, prob in t2_topics])\
    ax[1].set_title(text2_name)\
    \
    plt.tight_layout()\
    return fig\
\
def term_similarity_heatmap(text1, text2, max_terms=30):\
    # Get common important terms\
    vectorizer = TfidfVectorizer(max_features=max_terms)\
    X = vectorizer.fit_transform([text1, text2])\
    terms = vectorizer.get_feature_names()\
    \
    # Create term-term similarity matrix\
    similarity = X.T.dot(X).toarray()\
    \
    # Plot heatmap\
    plt.figure(figsize=(14, 12))\
    sns.heatmap(similarity, annot=True, fmt=".2f", \
               xticklabels=terms, yticklabels=terms)\
    plt.title("Term Similarity Between Annotations")\
    return plt.gcf()\
\
def plot_term_frequencies(text1, text2, top_n=30, text1_name = 'Freq1_', text2_name = 'Freq2_'):\
    # Get word frequencies\
    freq1 = Counter(text1.lower().split())\
    freq2 = Counter(text2.lower().split())\
    \
    # Create comparison dataframe\
    terms = set(list(freq1.keys())[:top_n] + list(freq2.keys())[:top_n])\
    df = pd.DataFrame(\{\
        'Term': list(terms),\
        text1_name: [freq1.get(term, 0) for term in terms],\
        text2_name: [freq2.get(term, 0) for term in terms]\
    \})\
    \
    # Plot\
    plt.figure(figsize=(14, 5))\
    sns.barplot(data=df.melt(id_vars='Term'), x='Term', y='value', hue='variable')\
    plt.xticks(rotation=90)\
    plt.tight_layout()\
\
def create_term_network(text1, text2, min_freq=1, max_terms=100):\
    # Count term frequencies\
    terms1 = Counter(text1.lower().split())\
    terms2 = Counter(text2.lower().split())\
    \
    # Combine and filter terms\
    all_terms = Counter()\
    all_terms.update(terms1)\
    all_terms.update(terms2)\
    \
    # Filter to most common terms that appear at least min_freq times\
    filtered_terms = [term for term, count in all_terms.most_common(max_terms) \
                     if count >= min_freq]\
    \
    # Create graph\
    G = nx.Graph()\
    \
    # Add nodes with attributes\
    for term in filtered_terms:\
        G.add_node(term, \
                  count1=terms1.get(term, 0),\
                  count2=terms2.get(term, 0),\
                  total=all_terms.get(term, 0))\
    \
    # Add edges based on co-occurrence\
    window_size = 5\
    \
    # Process text1\
    text1_words = text1.lower().split()\
    for i in range(len(text1_words)):\
        term1 = text1_words[i]\
        if term1 not in filtered_terms:\
            continue\
            \
        # Look at nearby terms within window\
        for j in range(max(0, i-window_size), min(len(text1_words), i+window_size+1)):\
            if i == j:\
                continue\
                \
            term2 = text1_words[j]\
            if term2 in filtered_terms:\
                if G.has_edge(term1, term2):\
                    G[term1][term2]['weight'] = G[term1][term2].get('weight', 0) + 1\
                else:\
                    G.add_edge(term1, term2, weight=1)\
    \
    # Repeat for text2\
    text2_words = text2.lower().split()\
    for i in range(len(text2_words)):\
        term1 = text2_words[i]\
        if term1 not in filtered_terms:\
            continue\
            \
        for j in range(max(0, i-window_size), min(len(text2_words), i+window_size+1)):\
            if i == j:\
                continue\
                \
            term2 = text2_words[j]\
            if term2 in filtered_terms:\
                if G.has_edge(term1, term2):\
                    G[term1][term2]['weight'] = G[term1][term2].get('weight', 0) + 1\
                else:\
                    G.add_edge(term1, term2, weight=1)\
    \
    # Remove edges with low weights\
    edges_to_remove = [(u, v) for u, v, d in G.edges(data=True) if d['weight'] < 2]\
    G.remove_edges_from(edges_to_remove)\
    \
    return G\
\
def plot_term_network(G, title="Term Network"):\
    # Get node attributes for sizing and coloring\
    # node_size = [G.nodes[n].get('total', 1) * 10 for n in G.nodes()]\
    node_size = [G.nodes[n].get('total', 1) for n in G.nodes()]\
    \
    # Calculate color based on which text it appears more in\
    node_color = []\
    for n in G.nodes():\
        count1 = G.nodes[n].get('count1', 0)\
        count2 = G.nodes[n].get('count2', 0)\
        \
        if count1 > count2:\
            # More in text1 - blue\
            color = (0, 0, 1, count1/(count1+count2))\
        elif count2 > count1:\
            # More in text2 - red\
            color = (1, 0, 0, count2/(count1+count2))\
        else:\
            # Equal - purple\
            color = (0.5, 0, 0.5, 0.7)\
            \
        node_color.append(color)\
    \
    # Plot\
    plt.figure(figsize=(12, 12))\
    pos = nx.spring_layout(G, k=0.9)\
    \
    max_weight = max([G[u][v].get('weight', 1) for u, v in G.edges()]) if G.edges() else 1\
    edge_colors = [(0.5, 0.5, 0.5, min(0.8, G[u][v].get('weight', 1)/max_weight)) for u, v in G.edges()]\
    \
    nx.draw_networkx(\
        G, pos=pos,\
        node_size=node_size,\
        node_color=node_color,\
        with_labels=True,\
        width=[G[u][v].get('weight', 1)/len(node_size) for u, v in G.edges()],\
        alpha=0.9,\
        # edge_color=edge_colors,\
        font_size=10\
    )\
    \
    plt.title(title)\
    plt.axis('off')\
    return plt.gcf()\
\
def plot_term_overlap(text1, text2):\
    terms1 = set(text1.lower().split())\
    terms2 = set(text2.lower().split())\
    \
    venn2([terms1, terms2], set_labels=('Comparison 1', 'Comparison 2'))\
\
def create_comparative_wordcloud(text1, text2, title1="Comparison 1", title2="Comparison 2"):\
    # Create word clouds for each text\
    wc1 = WordCloud(width=800, height=400, background_color='white').generate(text1)\
    wc2 = WordCloud(width=800, height=400, background_color='white').generate(text2)\
    \
    # Create a comparison plot\
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20,10))\
    ax1.imshow(wc1)\
    ax1.set_title(title1)\
    ax1.axis('off')\
    ax2.imshow(wc2)\
    ax2.set_title(title2)\
    ax2.axis('off')\
    plt.tight_layout()\
    return fig\
\
# For finding unique/shared terms\
def compare_terms(text1, text2):\
    # Convert texts to term sets\
    terms1 = set(text1.lower().split())\
    terms2 = set(text2.lower().split())\
    \
    unique_to_1 = terms1 - terms2\
    unique_to_2 = terms2 - terms1\
    shared = terms1 & terms2\
    \
    return unique_to_1, unique_to_2, shared\
\
\
def create_wordcloud(text1, title=""):\
    # Create word clouds for each text\
    wc1 = WordCloud(width=800, height=400, background_color='white').generate(text1)\
    \
    # Create a comparison plot\
    fig, ax1 = plt.subplots(1, 1, figsize=(20,10))\
    ax1.imshow(wc1)\
    ax1.set_title(title)\
    ax1.axis('off')\
    plt.tight_layout()\
    return fig\
\
def quick_split(term_str):\
    return term_str.split(' ')\
\
def unique_split(term_str):\
    return(' '.join(list(set(term_str.split(' ')))))\
\
additional_stop_words = \{*[str(x) + '.' for x in range(20)], \
                         *[str(x) + ')' for x in range(20)],\
                         *[str(x) + ':' for x in range(20)],\
                         *[str(x) for x in range(20)],\
                         'pathway', 'gene', 'pathways', 'genes', 'protein', 'proteins', 'cell', '**', '-', '.', ':', '(', ')', ',',\
                         'response', 'responses', 'role', 'involved', 'cellular', 'metabolism', 'metabolic','signaling','regulation','processing','regulates',\
                        'suggest','maintaining','processes','potentially','family','signal','related','potential','relevant','associated','regulators','regulate',\
                        'linked','components','identified','identify','identifies','influencing','influence','like','member','presence','functions','set','key','regulatory',\
                        'involvement','type','small','known','processing','involvement','plays','play','possibly','implicated','particularly','biology','converting',\
                        'suggests','movement','free','regulation','indicating','indicate','wound','players','balance','roles','coordinated','crucial','converts','energy',\
                        'receptor','retaining','activity','critical','controlling','various','residues','additionally','activity','directly','target','responsible','conversion',\
                        'responsible','expression','steps','binding','factor','containing','regulation','regulating','processing','including','modulating','cells','targeted','therapy',\
                        'regulation','processing','regulation','notably','impact','interplay','converge','diverse','group','points','properties','impact','biological','molecule',\
                        'channel','facilitate','interacts','volume','response','emphasis','control','modifies','regulation','modulation','interconnected','facilitates','indicates',\
                        'suggesting','process','organization','affecting','implications','modulate','integrity','survival','factors','supports','nature','important','flexibility',\
                        'tissues','strength','provides','playing','long','major','additional','shape','participate','participates','contribute','number','primarily','component',\
                        'based','function','structure','factors','collectively','sorting','contributing','changes','support','providing','supporting','contribute','initiation','wide',\
                        'negative','influences','progression','force','provides','maintains','ensuring','handling','generate','overall','heavy','highlights','indirectly','growth',\
                        'necessary','molecules','active','element','core','function','stress','specific','direct','generally','subsequent','enhancing','complex','range','production',\
                        'maintenance','miscellaneous','face','effort','detection','induced','presentation','attracting','focuses','integral','proper','conditions','link','effect',\
                        'linking','relationship','cycle','summary','responds','repeat','tissue','underscores','significant','adaptation','ensures','integrated','focus','mediating',\
                        'maintain','sites','affect','generates','interacting','unite','solute','modulates','interactions','rna','mrna'\}\
\
stop_words = set(additional_stop_words).union(set(ENGLISH_STOP_WORDS))\
\
def preprocess_terms(ora_res, unique = False):\
    if unique:\
        cleaned = ' '.join(list(set(re.split(r"[ |**|\\.|(|)|,]", ora_res.lower())) - set(stop_words)))\
    else:\
        temp = list(re.split(r"[ |**|\\.|(|)|,]", ' '.join(ora_res.lower().split('name'))))\
        temp_filt = [x for x in temp if x not in list(stop_words)]\
        cleaned = ' '.join(temp_filt)\
    return cleaned\
\
%who_ls function\
\
import numpy as np\
import matplotlib.pyplot as plt\
from wordcloud import WordCloud\
from collections import Counter\
import re\
\
def create_comparative_wordcloud_color(text1, text2, source1_name="Source 1", source2_name="Source 2", \
                                 color1='#ff9999', color2='#87cefa', common_color='#5e3ea5',\
                                 width=800, height=400, max_words=200, \
                                 background_color='white', min_font_size=8):\
    """\
    Create a color-coded wordcloud from two text sources.\
    \
    Parameters:\
    -----------\
    text1, text2 : str\
        The text content from each source\
    source1_name, source2_name : str\
        Names to label each source in the legend\
    color1, color2 : str\
        Colors for words unique to each source\
    common_color : str\
        Color for words that appear in both sources\
    width, height : int\
        Dimensions of the wordcloud\
    max_words : int\
        Maximum number of words to include\
    background_color : str\
        Background color of the wordcloud\
    min_font_size : int\
        Minimum font size for words\
        \
    Returns:\
    --------\
    fig : matplotlib.figure.Figure\
        The figure containing the wordcloud and legend\
    """\
    \
    # Process texts\
    def clean_text(text):\
        # Convert to lowercase and remove punctuation\
        text = re.sub(r'[^\\w\\s]', '', text.lower())\
        # Split into words\
        return text.split()\
    \
    # Clean and tokenize\
    words1 = clean_text(text1)\
    words2 = clean_text(text2)\
    \
    # Count word frequencies\
    word_counts1 = Counter(words1)\
    word_counts2 = Counter(words2)\
    \
    # Identify words unique to each source and common words\
    unique_words1 = set(word_counts1.keys()) - set(word_counts2.keys())\
    unique_words2 = set(word_counts2.keys()) - set(word_counts1.keys())\
    common_words = set(word_counts1.keys()) & set(word_counts2.keys())\
    \
    # Combine word frequencies\
    all_words = \{\}\
    \
    # Function to determine color based on word source\
    def color_func(word, font_size, position, orientation, random_state=None, **kwargs):\
        if word in unique_words1:\
            return color1\
        elif word in unique_words2:\
            return color2\
        else:\
            return common_color\
    \
    # Add all words with their frequencies\
    for word in word_counts1:\
        all_words[word] = word_counts1[word]\
    \
    for word in word_counts2:\
        if word in all_words:\
            # For common words, take the maximum frequency to give them proper weight\
            all_words[word] = max(all_words[word], word_counts2[word])\
        else:\
            all_words[word] = word_counts2[word]\
    \
    # Create wordcloud\
    wc = WordCloud(width=width, height=height, \
                   background_color=background_color,\
                   max_words=max_words,\
                   min_font_size=min_font_size,\
                   color_func=color_func,\
                   random_state=42).generate_from_frequencies(all_words)\
    \
    # Create plot\
    fig, ax = plt.subplots(figsize=(width/100, height/100))\
    ax.imshow(wc, interpolation='bilinear')\
    ax.axis('off')\
    \
    # Add legend with colored patches\
    import matplotlib.patches as mpatches\
    \
    #legend_patches = [\
    #    mpatches.Patch(color=color1, label=f'Unique to \{source1_name\}'),\
    #    mpatches.Patch(color=color2, label=f'Unique to \{source2_name\}'),\
    #    mpatches.Patch(color=common_color, label='Common to both')\
    #]\
    \
    # Place legend in the upper right corner\
    #plt.legend(handles=legend_patches, loc='upper right', \
    #           frameon=True, framealpha=0.8, facecolor='white')\
    \
    plt.title(f'Comparative Word Cloud: \{source1_name\} vs \{source2_name\}')\
    plt.tight_layout(pad=0)\
    \
    return fig\
\
L_ICM_fib_process = 'Extracellular matrix organization and immune response modulation\\\
1. **Extracellular Matrix (ECM) Components and Remodeling**: FN1 (fibronectin), VCAN (versican), COL1A2 (collagen type I alpha 2), and ASPN (asporin) are key components of the extracellular matrix. These proteins are involved in the structural organization and remodeling of the ECM, which is crucial for tissue integrity and repair. FN1 and VCAN play roles in cell adhesion, migration, and proliferation, while COL1A2 provides tensile strength to tissues. ASPN modulates collagen fibrillogenesis and can influence ECM mineralization.\\\
2. **Immune Response Modulation**: IL32 is a cytokine involved in the regulation of immune responses, including the activation of natural killer cells and the production of pro-inflammatory cytokines. SERPINE2, a serine protease inhibitor, can modulate inflammation and tissue remodeling by inhibiting proteases involved in ECM degradation. These proteins collectively contribute to the modulation of immune responses and inflammation, particularly in the context of tissue injury and repair.\\\
3. **Signal Transduction and Transcriptional Regulation**: NR4A1 and NR4A2 are nuclear receptors that function as transcription factors involved in the regulation of gene expression in response to various stimuli, including stress and inflammation. STAT3 is a transcription factor activated by cytokines and growth factors, playing a critical role in cell growth, survival, and immune responses. These proteins integrate extracellular signals to modulate cellular responses, including those related to ECM organization and immune function.\\\
4. **Protein Degradation and Turnover**: PSMD14 and PSMD2 are components of the 26S proteasome, involved in the degradation of ubiquitinated proteins. This process is essential for maintaining protein homeostasis and regulating the turnover of proteins involved in ECM remodeling and immune responses.\\\
5. **Cytoskeletal Dynamics and Cell Motility**: TUBA1B (tubulin alpha-1B chain) and PDLIM3 are involved in cytoskeletal organization and cell motility. These proteins contribute to cellular processes such as migration and adhesion, which are important for tissue repair and immune cell function.\\\
Overall, the system of interacting proteins is primarily involved in the organization and remodeling of the extracellular matrix, as well as the modulation of immune responses, particularly in the context of tissue repair and inflammation. The interplay between ECM components, immune modulators, and signaling pathways highlights the integrated nature of these biological processes.'\
\
L_ICM_fib_bio = 'The proteins identified can be grouped into several functional categories, suggesting potential biological processes and disease associations:\\\
1. **Transcription Factors and Nuclear Receptors**: NR4A1, NR4A2, and NFIL3 are nuclear receptors and transcription factors involved in regulating gene expression. NR4A1 and NR4A2 are part of the NR4A subfamily, which is implicated in cellular stress responses, apoptosis, and metabolism. NFIL3 is involved in circadian rhythm regulation and immune responses. These proteins may collectively influence gene expression patterns in response to environmental or cellular signals.\\\
2. **Extracellular Matrix and Cell Adhesion**: FN1 (fibronectin), VCAN (versican), COL1A2 (collagen type I alpha 2), and ASPN (asporin) are components of the extracellular matrix (ECM) and are involved in cell adhesion, migration, and tissue remodeling. SERPINE2, a serine protease inhibitor, regulates ECM degradation. These proteins suggest a role in tissue remodeling, wound healing, or cancer metastasis.\\\
3. **Immune Response and Inflammation**: IL32 is a cytokine involved in inflammatory responses. B2M (beta-2-microglobulin) and PSMB8 (proteasome subunit beta type-8) are involved in antigen processing and presentation. These proteins indicate a potential role in immune regulation and inflammatory diseases.\\\
4. **Signal Transduction and Cell Cycle Regulation**: STAT3 is a transcription factor activated by cytokines and growth factors, playing a role in cell growth and apoptosis. CDK2AP1 is involved in cell cycle regulation. These proteins may be part of pathways controlling cell proliferation and survival, possibly linked to cancer.\\\
5. **Protein Synthesis and Degradation**: RPL36A and RPL17 are ribosomal proteins involved in protein synthesis. PSMD14 and PSMD2 are components of the proteasome, involved in protein degradation. These proteins suggest a role in maintaining protein homeostasis.\\\
6. **Cytoskeletal and Structural Proteins**: TUBA1B (tubulin alpha-1B) and PDLIM3 are involved in cytoskeletal organization and muscle function. These proteins may play a role in maintaining cellular structure and motility.\\\
7. **Metabolic and Mitochondrial Proteins**: NDUFA13 is part of the mitochondrial electron transport chain, indicating a role in cellular energy metabolism. MGST3 is involved in detoxification processes.\\\
8. **RNA Processing and Regulation**: DCP2 is involved in mRNA decapping, a critical step in mRNA degradation. XIST is a long non-coding RNA involved in X-chromosome inactivation. These proteins suggest roles in RNA stability and gene expression regulation.\\\
The identified proteins suggest involvement in processes such as ECM remodeling, immune response, cell cycle regulation, and protein homeostasis. The presence of ECM components and immune-related proteins indicates potential roles in tissue remodeling and inflammatory diseases. The involvement of transcription factors and cell cycle regulators suggests possible links to cancer. The proteins related to RNA processing and mitochondrial function highlight roles in gene expression regulation and energy metabolism. Overall, these proteins may collectively participate in complex biological processes, including cancer progression, immune regulation, and tissue remodeling.'\
\
L_AMI_fib_process = 'Extracellular matrix organization and ribosomal function\\\
1. **Extracellular Matrix (ECM) Organization**: COL1A2, BGN, SPARC, AEBP1, and VCAN are involved in the structure and function of the extracellular matrix. COL1A2 encodes a component of type I collagen, a major structural protein in the ECM. BGN (biglycan) and VCAN (versican) are proteoglycans that interact with collagen and other ECM components, influencing tissue resilience and cell signaling. SPARC (secreted protein acidic and rich in cysteine) modulates cell-matrix interactions and ECM assembly. AEBP1 (adipocyte enhancer-binding protein 1) is implicated in ECM remodeling and regulation of adipogenesis.\\\
2. **Ribosomal Function and Protein Synthesis**: RPS4Y1, RPS26, and RPS4X are ribosomal proteins that are integral to the structure and function of ribosomes, facilitating protein synthesis. These proteins are components of the small ribosomal subunit and play roles in mRNA translation and ribosome assembly.\\\
3. **RNA Processing and Regulation**: DDX3X and DDX3Y are RNA helicases involved in RNA processing, translation initiation, and regulation of gene expression. They participate in unwinding RNA structures, which is crucial for various aspects of RNA metabolism.\\\
4. **Cellular Signaling and Regulation**: LATS2 is a serine/threonine kinase involved in the Hippo signaling pathway, which regulates cell proliferation and apoptosis. CD99 is involved in cell adhesion and migration, playing roles in immune response and cancer progression.\\\
5. **Cytoskeletal Dynamics**: TUBA1B encodes a tubulin protein that is a key component of microtubules, essential for maintaining cell shape, intracellular transport, and cell division.\\\
6. **Miscellaneous Functions**: XIST is a long non-coding RNA involved in X-chromosome inactivation, a process crucial for dosage compensation in females. HMGN3 is a chromatin-associated protein that modulates chromatin structure and gene expression.\\\
The proteins in this system are primarily involved in two distinct biological processes: ECM organization and ribosomal function. The presence of proteins involved in RNA processing, cellular signaling, and cytoskeletal dynamics suggests additional roles in cellular regulation and structure. However, the ECM organization and ribosomal function are the most prominent processes based on the number of proteins involved.'\
\
L_AMI_fib_bio = '1. **Ribosomal Proteins and Translation Regulation**:\\\
   - **RPS4Y1, RPS26, RPS4X**: These proteins are components of the ribosome, involved in protein synthesis. RPS4Y1 and RPS4X are homologous, with RPS4Y1 being Y-linked and RPS4X being X-linked, suggesting a role in sex-specific differences in translation. RPS26 is a component of the small ribosomal subunit, essential for mRNA translation.\\\
   - **DDX3Y, DDX3X**: These are DEAD-box RNA helicases involved in RNA metabolism, including translation initiation. DDX3X is X-linked, while DDX3Y is Y-linked, and both may have roles in sex-specific gene expression regulation.\\\
2. **Extracellular Matrix and Cell Adhesion**:\\\
   - **COL1A2, BGN, AEBP1, SPARC, VCAN**: These proteins are involved in the structure and function of the extracellular matrix (ECM). COL1A2 is a collagen component, BGN (biglycan) and VCAN (versican) are proteoglycans, SPARC (osteonectin) is involved in ECM remodeling, and AEBP1 is associated with ECM organization. These proteins collectively contribute to tissue integrity and cell adhesion, potentially influencing processes like wound healing and fibrosis.\\\
3. **Cell Cycle and Growth Regulation**:\\\
   - **LATS2**: A key component of the Hippo signaling pathway, which regulates cell proliferation and apoptosis. It might interact with other proteins to control cell growth and prevent tumorigenesis.\\\
   - **CD99**: Involved in cell adhesion and apoptosis, potentially interacting with LATS2 in growth regulation contexts.\\\
4. **Protease Activity and Inhibition**:\\\
   - **PRSS23, SERPINE2**: PRSS23 is a serine protease, while SERPINE2 is a serine protease inhibitor. Their interaction suggests a regulatory mechanism in proteolytic processes, which could be relevant in tissue remodeling or inflammation.\\\
5. **Microtubule Dynamics**:\\\
   - **TUBA1B**: A tubulin protein involved in microtubule formation, essential for cell division and intracellular transport.\\\
6. **Non-coding RNA and Gene Regulation**:\\\
   - **XIST**: A long non-coding RNA crucial for X-chromosome inactivation in females, indicating a role in dosage compensation and gene expression regulation.\\\
7. **Chromatin and Gene Expression**:\\\
   - **HMGN3**: A chromatin-associated protein that modulates chromatin structure and gene expression, potentially interacting with other regulatory proteins to influence transcription.\\\
8. **Pancreatic Development and Function**:\\\
   - **PPDPF**: Involved in pancreatic development, though its specific interactions in this context are less clear.\\\
The identified proteins suggest involvement in several key biological processes, including translation regulation, ECM organization, cell cycle control, and protease activity. The presence of ribosomal proteins and RNA helicases indicates a focus on protein synthesis and RNA metabolism, with potential sex-specific regulatory roles due to the X and Y-linked homologs. ECM-related proteins highlight roles in tissue structure and repair, while LATS2 and CD99 suggest regulation of cell growth and apoptosis. The interplay between proteases and inhibitors points to controlled proteolysis, relevant in tissue remodeling. XIST and HMGN3 indicate gene expression regulation, with XIST specifically involved in X-chromosome inactivation. While no single disease process is immediately apparent, the combination of these proteins could be relevant in conditions involving aberrant cell growth, ECM dysregulation, or sex-specific gene expression differences.'\
\
L_ICM_mac_process = 'Immune Response Modulation and Inflammatory Regulation  \\\
1. **CD74** acts as a receptor for macrophage migration inhibitory factor (MIF) and is involved in antigen presentation, playing a crucial role in immune response regulation. **MIF** itself is a pro-inflammatory cytokine that modulates immune responses and is involved in the regulation of macrophage function.\\\
2. **TGFBI** (Transforming Growth Factor Beta Induced) is involved in cell adhesion and migration, often associated with immune response modulation and tissue repair processes. **KLF4** and **KLF6** are transcription factors that regulate inflammation and immune responses, with KLF4 being involved in macrophage polarization.\\\
3. **PSAP** (Prosaposin) is involved in lysosomal function and has roles in modulating immune responses. **GABARAP** is associated with autophagy, a process that can influence immune responses by degrading intracellular pathogens.\\\
4. **MSR1** (Macrophage Scavenger Receptor 1) is involved in the recognition and clearance of pathogens, playing a role in innate immunity. **FPR3** (Formyl Peptide Receptor 3) is involved in chemotaxis and immune cell signaling.\\\
5. **CCL18** is a chemokine involved in attracting immune cells to sites of inflammation, while **CD14** is a co-receptor for the detection of bacterial lipopolysaccharide, playing a role in the innate immune response.\\\
6. **ADAMTS2** is involved in extracellular matrix remodeling, which can influence immune cell migration and tissue repair. **APOD** (Apolipoprotein D) has roles in lipid metabolism and may influence inflammatory processes.\\\
7. **SLC40A1** (Ferroportin) is involved in iron homeostasis, which can affect immune cell function and inflammatory responses. **RGL1** and **CFD** (Complement Factor D) are involved in signaling pathways that can modulate immune responses.\\\
8. **AHNAK** is involved in cellular signaling and may influence immune cell function. **PLCG2** (Phospholipase C Gamma 2) is involved in signal transduction pathways that regulate immune cell activation.\\\
The proteins in this system are primarily involved in modulating immune responses and regulating inflammation, with roles in antigen presentation, cytokine signaling, chemotaxis, and immune cell activation. The interplay between these proteins suggests a coordinated effort to regulate immune responses and maintain homeostasis in the face of inflammatory stimuli.'\
\
L_ICM_mac_bio = 'The set of proteins provided can be grouped into several functional categories based on their known roles and interactions:\\\
1. **Immune Response and Inflammation**: \\\
   - **CD74, MSR1, FPR3, CCL18, CD14, LILRB5, MIF**: These proteins are involved in immune response and inflammation. CD74 acts as a receptor for macrophage migration inhibitory factor (MIF) and is involved in antigen presentation. MSR1 is a scavenger receptor involved in the immune response. FPR3 is a receptor involved in chemotaxis of immune cells. CCL18 is a chemokine involved in immune cell recruitment. CD14 is a co-receptor for the detection of bacterial lipopolysaccharide. LILRB5 is an inhibitory receptor on immune cells. These proteins suggest a coordinated role in modulating immune responses and inflammation.\\\
2. **Metabolism and Energy Regulation**:\\\
   - **PDK4, SLC40A1, APOD, GATM**: PDK4 is involved in the regulation of glucose metabolism by inhibiting the pyruvate dehydrogenase complex. SLC40A1 is involved in iron transport. APOD is associated with lipid metabolism. GATM is involved in creatine biosynthesis. These proteins suggest a role in metabolic regulation, possibly in response to cellular stress or energy demands.\\\
3. **Cellular Stress Response and Apoptosis**:\\\
   - **TSC22D3, PSAP, AHNAK, MTRNR2L12, MTRNR2L8, RBM3**: TSC22D3 is a glucocorticoid-induced leucine zipper protein involved in stress response. PSAP is involved in lysosomal function and apoptosis. AHNAK is a large scaffold protein involved in cellular architecture and stress response. MTRNR2L12 and MTRNR2L8 are mitochondrial-derived peptides with potential roles in apoptosis. RBM3 is a cold-shock protein involved in stress response. These proteins may be involved in cellular adaptation to stress and regulation of apoptosis.\\\
4. **Signal Transduction and Cell Communication**:\\\
   - **PLCG2, IQGAP2, KLF4, KLF6, RGL1**: PLCG2 is involved in signal transduction pathways, particularly in immune cells. IQGAP2 is a scaffold protein involved in cell signaling and cytoskeletal organization. KLF4 and KLF6 are transcription factors involved in cell differentiation and proliferation. RGL1 is involved in Ras signaling pathways. These proteins suggest roles in modulating cell signaling and communication.\\\
5. **Protein Synthesis and Ribosomal Function**:\\\
   - **RPL36A, RPS10, MRPL23, EIF1AY**: These proteins are components of the ribosome, involved in protein synthesis. RPL36A and RPS10 are part of the cytoplasmic ribosome, while MRPL23 is part of the mitochondrial ribosome. EIF1AY is involved in the initiation of translation. These proteins indicate a role in protein synthesis and possibly in the regulation of translation under specific conditions.'\
\
L_AMI_mac_process = 'Inflammatory Response and Metabolic Regulation\\\
1. **HIF1A (Hypoxia-Inducible Factor 1-alpha)** is a transcription factor that responds to changes in available oxygen in the cellular environment, playing a crucial role in cellular adaptation to hypoxia. It regulates genes involved in energy metabolism, angiogenesis, and apoptosis, which are essential for inflammatory responses and metabolic regulation.\\\
2. **CXCL8 (Interleukin-8)** is a chemokine involved in the inflammatory response, primarily acting as a chemoattractant for neutrophils. It plays a significant role in the recruitment and activation of immune cells during inflammation.\\\
3. **CD36** is a scavenger receptor involved in fatty acid metabolism and inflammation. It facilitates the uptake of long-chain fatty acids and oxidized low-density lipoproteins, linking lipid metabolism with inflammatory processes.\\\
4. **CTSB (Cathepsin B)** is a lysosomal cysteine protease involved in protein degradation and processing. It participates in the regulation of inflammation and apoptosis, contributing to tissue remodeling and immune responses.\\\
5. **HMOX1 (Heme Oxygenase 1)** is an enzyme that degrades heme into biliverdin, free iron, and carbon monoxide. It has anti-inflammatory and cytoprotective effects, modulating oxidative stress and inflammation.\\\
6. **FOXO3** is a transcription factor involved in the regulation of oxidative stress, apoptosis, and metabolism. It plays a role in modulating inflammatory responses and maintaining cellular homeostasis.\\\
7. **FABP5 (Fatty Acid Binding Protein 5)** is involved in the intracellular transport of fatty acids and regulation of lipid metabolism. It is implicated in inflammatory processes and energy homeostasis.\\\
8. **ADIPOR2 (Adiponectin Receptor 2)** is a receptor for adiponectin, a hormone involved in glucose regulation and fatty acid oxidation. It plays a role in anti-inflammatory signaling and metabolic regulation.\\\
The proteins in this system are primarily involved in the regulation of inflammatory responses and metabolic processes. The interplay between these proteins highlights the integration of metabolic regulation with immune responses, where factors like HIF1A and FOXO3 modulate cellular responses to stress, while CXCL8 and CD36 link inflammation with metabolic pathways. This system underscores the complex relationship between inflammation and metabolism, with potential implications for conditions such as metabolic syndrome and chronic inflammatory diseases.'\
\
L_AMI_mac_bio = 'Hypoxia and Inflammatory Response\\\
- **HIF1A (Hypoxia-Inducible Factor 1-alpha)**: A key regulator of the cellular response to low oxygen levels, influencing genes involved in energy metabolism, angiogenesis, and apoptosis.\\\
- **CXCL8 (Interleukin-8)**: A chemokine involved in the inflammatory response, attracting neutrophils to sites of infection or injury.\\\
- **HMOX1 (Heme Oxygenase 1)**: Provides cytoprotection against oxidative stress and is upregulated by hypoxia and inflammation.\\\
- **MAP3K8 (Mitogen-Activated Protein Kinase Kinase Kinase 8)**: Involved in the MAPK signaling pathway, which is activated by inflammatory cytokines and stress.\\\
These proteins suggest a coordinated response to hypoxic and inflammatory conditions. HIF1A can upregulate CXCL8 and HMOX1, promoting inflammation and protection against oxidative stress. MAP3K8 may further propagate inflammatory signaling.\\\
 Lipid Metabolism and Transport\\\
- **CD36**: A receptor involved in fatty acid uptake and lipid metabolism.\\\
- **FABP5 (Fatty Acid Binding Protein 5)**: Facilitates intracellular transport of fatty acids.\\\
- **ADIPOR2 (Adiponectin Receptor 2)**: Mediates the effects of adiponectin, enhancing fatty acid oxidation and glucose uptake.\\\
These proteins are involved in lipid metabolism, suggesting a role in energy homeostasis. CD36 and FABP5 work together to manage fatty acid transport and utilization, while ADIPOR2 modulates metabolic pathways.\\\
Protein Degradation and Ubiquitination\\\
- **KLHL18**: A substrate adaptor for a cullin-RING E3 ubiquitin ligase complex, involved in protein ubiquitination.\\\
- **RNF44 (Ring Finger Protein 44)**: An E3 ubiquitin-protein ligase, contributing to protein degradation.\\\
- **SMARCC1**: Part of the SWI/SNF chromatin remodeling complex, indirectly influencing protein degradation pathways.\\\
These proteins are involved in protein ubiquitination and degradation, suggesting a role in maintaining protein homeostasis and regulating protein levels in response to cellular signals.\\\
Cell Cycle and Apoptosis\\\
- **FOXO3 (Forkhead Box O3)**: A transcription factor involved in cell cycle regulation and apoptosis.\\\
- **GRAMD4**: A mitochondrial protein that can induce apoptosis.\\\
- **RHOB**: A small GTPase involved in apoptosis and cytoskeletal organization.\\\
These proteins suggest a role in regulating cell survival and apoptosis. FOXO3 and GRAMD4 may promote apoptosis under stress conditions, while RHOB can influence cytoskeletal dynamics and cell death.\\\
The identified proteins suggest involvement in hypoxia and inflammatory responses, lipid metabolism, protein degradation, and cell cycle regulation. The interplay between HIF1A, CXCL8, and HMOX1 indicates a response to hypoxic and inflammatory conditions, potentially relevant to diseases like cancer and chronic inflammatory disorders. The lipid metabolism group, including CD36 and FABP5, suggests roles in metabolic diseases such as obesity and diabetes. The protein degradation group highlights mechanisms for maintaining protein homeostasis, while the cell cycle and apoptosis group points to pathways regulating cell survival, potentially relevant to cancer and neurodegenerative diseases.'\
\
K_IZ_fib_process = 'Name: Cardiac Muscle Contraction and Structural Integrity\\\
1. **Contractile Proteins**: The proteins TNNT2, ACTA1, TTN, DES, MYL3, TNNC1, MYL2, MYH7, MYH6, and ACTG1 are integral components of the sarcomere, the fundamental unit of muscle contraction. TNNT2 (troponin T) and TNNC1 (troponin C) are part of the troponin complex that regulates muscle contraction in response to calcium ions. MYH7 and MYH6 are myosin heavy chains that interact with actin filaments (ACTA1 and ACTG1) to generate force. TTN (titin) provides structural support and elasticity, while DES (desmin) maintains the structural integrity of the sarcomere.\\\
2. **Structural and Regulatory Proteins**: ANKRD1, CRYAB, HSPB1, and HSPB7 are involved in maintaining the structural integrity and stress response of cardiac muscle cells. ANKRD1 is a cardiac ankyrin repeat protein that interacts with titin and is involved in mechanotransduction. CRYAB (alpha B-crystallin) and HSPB1 (heat shock protein beta-1) function as molecular chaperones, protecting cardiac cells from stress-induced damage.\\\
3. **Calcium Handling and Signaling**: Proteins such as MYLK (myosin light chain kinase) and TECRL (trans-2,3-enoyl-CoA reductase-like) are involved in calcium signaling pathways that regulate muscle contraction. MYLK phosphorylates myosin light chains, enhancing the interaction between actin and myosin.\\\
4. **Extracellular Matrix and Cell Adhesion**: FLNA (filamin A), FLNC (filamin C), and TAGLN (transgelin) contribute to the structural framework and mechanical stability of cardiac tissue. They interact with the cytoskeleton and extracellular matrix, facilitating cell adhesion and signal transduction.\\\
5. **Metabolic and Transport Proteins**: FABP3 (fatty acid-binding protein 3) and SLC27A3 (solute carrier family 27 member 3) are involved in fatty acid transport and metabolism, crucial for energy production in cardiac muscle cells.\\\
6. **Additional Proteins**: Proteins such as A2M (alpha-2-macroglobulin), AEBP1 (AE binding protein 1), and IGFBP4 (insulin-like growth factor-binding protein 4) may play roles in modulating extracellular matrix interactions and growth factor signaling, indirectly influencing cardiac muscle function.\\\
Overall, this system of interacting proteins is primarily involved in the contraction, structural integrity, and metabolic regulation of cardiac muscle, ensuring efficient heart function.'\
\
K_IZ_fib_bio = '## Group 1: S100 Proteins and Inflammation\\\
- **S100A8, S100A9, S100A4**: These proteins are part of the S100 family, known for their roles in inflammation and immune responses. S100A8 and S100A9 often form a heterodimer called calprotectin, which is involved in the regulation of inflammatory processes and immune cell chemotaxis. S100A4 is associated with metastasis and cancer progression.\\\
Muscle and Structural Proteins\\\
- **TNNT2, ACTA1, TTN, DES, MYL3, TNNC1, MYH7, MYH6, ACTG1, MYL2, MYL6, ACTA2, MYH11, FLNA, CNN1**: These proteins are primarily involved in muscle contraction and structural integrity. They are key components of the sarcomere, the fundamental unit of muscle contraction. Mutations in these proteins are often linked to cardiomyopathies and other muscle-related diseases.\\\
Heat Shock Proteins and Stress Response\\\
- **HSPA5, HSPB1, CRYAB, HSPB7**: These proteins are part of the heat shock protein family, which assists in protein folding and protection against stress-induced damage. They play roles in cellular stress responses and are implicated in neurodegenerative diseases and cancer.\\\
Immune and Inflammatory Mediators\\\
- **LYZ, CCL2, CCL21, MARCO, SRGN**: These proteins are involved in immune responses. LYZ (lysozyme) is an antimicrobial enzyme, while CCL2 and CCL21 are chemokines that recruit immune cells to sites of inflammation. MARCO is a scavenger receptor involved in pathogen recognition, and SRGN (serglycin) is involved in the storage and secretion of proteases in immune cells.\\\
Metabolic and Transport Proteins\\\
- **LAMTOR5, SLC27A3, SLC40A1, FABP3**: These proteins are involved in metabolic processes and transport. LAMTOR5 is part of the Ragulator complex, which regulates mTORC1 signaling. SLC27A3 and SLC40A1 are involved in fatty acid and iron transport, respectively. FABP3 is a fatty acid-binding protein important in lipid metabolism.\\\
Signaling and Regulatory Proteins\\\
- **NOTCH3, ECSCR, CALCRL, TIE1, HEG1, AOC3, ADAMTS1**: These proteins are involved in signaling pathways. NOTCH3 is part of the Notch signaling pathway, which is crucial for cell differentiation. ECSCR and CALCRL are involved in endothelial cell signaling. TIE1 and HEG1 are associated with vascular development. AOC3 is involved in amine oxidation and inflammation, while ADAMTS1 is a metalloprotease involved in extracellular matrix remodeling.\\\
Miscellaneous Proteins\\\
- **APP, GNL3, RNF181, MLXIP, DIPK2B, BCLAF1, MTRNR2L12, XIRP2, TECRL, MUSTN1, CLU, ADIRF, LBH, A2M, AEBP1, IGFBP4**: These proteins have diverse roles. APP is associated with Alzheimers disease. GNL3 and RNF181 are involved in cell cycle regulation. MLXIP is a transcription factor involved in glucose metabolism. BCLAF1 is involved in apoptosis regulation. CLU (clusterin) is implicated in neurodegenerative diseases and cancer.\\\
The identified proteins suggest a complex interplay between inflammation, muscle function, stress response, immune regulation, and metabolic processes. The presence of S100 proteins and chemokines indicates a strong inflammatory component, potentially linked to immune responses or inflammatory diseases. The muscle and structural proteins suggest involvement in muscle function or cardiomyopathies. Heat shock proteins indicate a response to cellular stress, which could be relevant in neurodegenerative diseases. The signaling proteins suggest roles in vascular development and cell differentiation. Overall, the proteins are associated with processes that could be relevant to inflammatory diseases, cardiomyopathies, neurodegenerative diseases, and cancer.'\
\
K_BZ_fib_process = 'Name: Cellular energy metabolism and regulation\\\
1. **PURA** (Purine-rich element binding protein A) is involved in the regulation of DNA replication and transcription. It plays a role in the control of cell growth and differentiation, and it has been implicated in the regulation of mitochondrial function, which is crucial for cellular energy metabolism.\\\
2. **CUL4A** (Cullin 4A) is a core component of the E3 ubiquitin-protein ligase complex, which is involved in the ubiquitination and subsequent proteasomal degradation of target proteins. This process is essential for maintaining cellular homeostasis and regulating various cellular processes, including cell cycle progression and DNA repair, which indirectly influence cellular energy metabolism.\\\
3. **SLC35F5** is a member of the solute carrier family, which is involved in the transport of molecules across cellular membranes. While its specific substrates and functions are not well-characterized, solute carriers generally play roles in cellular metabolism by facilitating the transport of metabolites and ions necessary for energy production and cellular function.\\\
4. **SORBS2** (Sorbin and SH3 domain-containing protein 2) is involved in cytoskeletal organization and signal transduction. It has been associated with insulin signaling pathways, which are critical for glucose uptake and energy metabolism in cells.\\\
5. **NDUFA6** is a subunit of the NADH:ubiquinone oxidoreductase (Complex I) in the mitochondrial electron transport chain. It plays a direct role in cellular energy production by participating in oxidative phosphorylation, a process that generates ATP, the primary energy currency of the cell.\\\
In summary, the proteins in this system are primarily involved in processes related to cellular energy metabolism and regulation. They contribute to mitochondrial function, energy production, and the regulation of metabolic pathways, highlighting their collective role in maintaining cellular energy homeostasis.'\
\
K_BZ_fib_bio = '**PURA and SORBS2**\\\
PURA (Purine-rich element binding protein A) is known for its role in DNA replication and transcription regulation. It is involved in neuronal development and has been implicated in neurodevelopmental disorders. SORBS2 (Sorbin and SH3 domain-containing protein 2) is involved in cytoskeletal organization and signal transduction, particularly in cardiac and skeletal muscle tissues. While these proteins have distinct primary functions, they both play roles in cellular signaling and structural organization. It might be the case that PURA and SORBS2 interact in a context involving cellular signaling pathways that affect both neuronal and muscular systems, potentially influencing developmental processes.\\\
**CUL4A and NDUFA6**\\\
CUL4A (Cullin 4A) is a core component of the E3 ubiquitin-protein ligase complex, which is involved in targeting proteins for degradation. It plays a role in DNA repair, cell cycle regulation, and transcriptional regulation. NDUFA6 (NADH:Ubiquinone Oxidoreductase Subunit A6) is a component of the mitochondrial respiratory chain complex I, essential for ATP production. While these proteins are involved in different cellular processes, they both contribute to cellular homeostasis and energy regulation. CUL4As role in protein degradation could influence mitochondrial function by regulating the turnover of mitochondrial proteins, including components of the respiratory chain like NDUFA6.\\\
**SLC35F5**\\\
SLC35F5 (Solute Carrier Family 35 Member F5) is a less characterized protein, but as a member of the solute carrier family, it is likely involved in the transport of molecules across cellular membranes. Its specific substrates and physiological roles are not well-defined, but it may play a role in cellular metabolism and homeostasis.\\\
The proteins identified\'97PURA, CUL4A, SLC35F5, SORBS2, and NDUFA6\'97suggest potential involvement in cellular signaling, structural organization, and metabolic regulation. PURA and SORBS2 may interact in pathways affecting neuronal and muscular development, while CUL4A and NDUFA6 could be linked through mechanisms of protein turnover and mitochondrial function. SLC35F5s role remains speculative, but it may contribute to cellular transport processes. There is no strong evidence directly linking these proteins to a specific disease process, but their functions suggest potential implications in neurodevelopmental disorders, muscular diseases, and metabolic conditions.'\
\
AF_process = 'Cytoskeletal Dynamics and Extracellular Matrix Remodeling\\\
1. **Cytoskeletal Dynamics**: Proteins such as MYH9, MYL9, MYL12A, MYL3, MYL2, TPM1, and TPM2 are involved in the regulation of actin-myosin interactions, which are crucial for muscle contraction and cellular motility. VIM and SPTAN1 contribute to the structural integrity and flexibility of the cytoskeleton, while CFL1 is involved in actin filament depolymerization, facilitating cytoskeletal reorganization.\\\
2. **Extracellular Matrix (ECM) Remodeling**: FN1, BGN, COL6A1, and COL6A2 are key components of the ECM, playing roles in cell adhesion, migration, and tissue integrity. TIMP1 regulates matrix metalloproteinases, which are involved in ECM degradation and remodeling. CCN2 is involved in ECM production and cell proliferation.\\\
3. **Signal Transduction and Regulation**: Proteins like PDGFRB and ARHGEF1 are involved in signal transduction pathways that regulate cell growth, migration, and cytoskeletal organization. RAB11B is involved in vesicular trafficking, which is essential for membrane dynamics and signal transduction.\\\
4. **Transcriptional Regulation and Stress Response**: JUNB, EGR1, and ID3 are transcription factors that regulate gene expression in response to various stimuli, including stress and growth signals. DUSP1 is involved in the dephosphorylation of MAP kinases, modulating cellular responses to stress.\\\
5. **Protein Synthesis and Ribosomal Function**: RPS9, RPL36, RPL41, RPL3L, and RPL35 are components of the ribosome, essential for protein synthesis. Their presence indicates active protein translation, which supports cellular growth and adaptation.\\\
The system of interacting proteins primarily highlights the interplay between cytoskeletal dynamics and ECM remodeling, with additional roles in signal transduction and transcriptional regulation. This integrated network supports cellular structure, motility, and response to environmental changes.'\
\
AF_bio = '### Protein Groups and Their Roles\\\
1. **Cytoskeletal and Muscle-Related Proteins:**\\\
   - **MYH9, MYL9, MYL12A, MYL3, MYL2, MYH7, TPM1, TPM2, TTN, MYBPC3, SPTAN1, VIM:** These proteins are involved in muscle contraction and cytoskeletal organization. MYH9, MYH7, and MYBPC3 are myosin heavy chains, while MYL9, MYL12A, MYL3, and MYL2 are myosin light chains, all crucial for muscle function. TPM1 and TPM2 are tropomyosins that regulate actin-myosin interactions. TTN (titin) is essential for muscle elasticity and structure. SPTAN1 (spectrin) and VIM (vimentin) are involved in maintaining cell integrity and cytoskeletal structure.\\\
   - **Potential Process:** These proteins suggest a role in muscle contraction and structural integrity, possibly indicating a focus on cardiac or skeletal muscle function.\\\
2. **Extracellular Matrix and Cell Adhesion:**\\\
   - **FN1, BGN, COL6A1, COL6A2, AEBP1, CCN2, TIMP1:** These proteins are involved in extracellular matrix (ECM) composition and remodeling. FN1 (fibronectin) and CCN2 (connective tissue growth factor) are key in cell adhesion and ECM interactions. BGN (biglycan) and TIMP1 (tissue inhibitor of metalloproteinases) regulate ECM stability and degradation. COL6A1 and COL6A2 are collagen components, crucial for ECM structure.\\\
   - **Potential Process:** This group suggests involvement in tissue repair, fibrosis, or ECM remodeling, which could be relevant in wound healing or fibrotic diseases.\\\
3. **Ribosomal and Translation-Related Proteins:**\\\
   - **RPS9, RPL36, RPL41, RPL3L, RPL35, RPL27, RPS8:** These ribosomal proteins are essential for protein synthesis. They are components of the ribosome, facilitating mRNA translation into proteins.\\\
   - **Potential Process:** This indicates active protein synthesis, possibly in response to increased cellular demands or stress.\\\
4. **Signal Transduction and Transcription Regulation:**\\\
   - **EGR1, JUNB, TCF25, TP53INP2, ID3, MAZ, TRIO, RASGRF1:** These proteins are involved in signal transduction and transcriptional regulation. EGR1 and JUNB are transcription factors that respond to growth signals. TP53INP2 is involved in p53-mediated apoptosis. TRIO and RASGRF1 are guanine nucleotide exchange factors, playing roles in signal transduction pathways.\\\
   - **Potential Process:** This group suggests regulation of cell growth, apoptosis, and response to external stimuli, potentially linking to cancer or stress responses.'\
\
CCR2_process = 'Extracellular Matrix Organization and Protein Synthesis\\\
1. **Extracellular Matrix (ECM) Components and Remodeling**: COL1A2, COL3A1, and COL6A2 are collagen proteins that are integral to the structure and function of the ECM, providing tensile strength and structural support. SPARC is involved in ECM remodeling and influences cell-matrix interactions. These proteins collectively contribute to the maintenance and organization of the ECM, which is crucial for tissue integrity and repair.\\\
2. **Ribosomal Proteins and Protein Synthesis**: RPL23, RPL13, RPLP1, and RPS26 are ribosomal proteins that play essential roles in the assembly and function of ribosomes, facilitating protein synthesis. TMA7 is involved in ribosome biogenesis, further supporting the process of translation.\\\
3. **Cellular Stress Response and Regulation**: DUSP1 is a dual-specificity phosphatase involved in the negative regulation of MAPK signaling pathways, which are activated in response to stress. TXNIP and SIRT6 are involved in cellular redox regulation and stress responses, with TXNIP acting as a regulator of oxidative stress and SIRT6 playing a role in DNA repair and metabolism.\\\
4. **Transport and Trafficking**: RAB21 and ARL8A are small GTPases involved in vesicular trafficking and endosomal transport, essential for maintaining cellular homeostasis and protein sorting.\\\
5. **Miscellaneous Functions**: A2M is a protease inhibitor involved in inhibiting a wide range of proteases, contributing to immune response and inflammation regulation. FABP5 is involved in fatty acid transport and metabolism, while SHANK3 is associated with synaptic function and neuronal signaling.\\\
The proteins in this system are primarily involved in ECM organization and protein synthesis, with additional roles in cellular stress response and intracellular transport. The interplay between these processes supports cellular structure, function, and adaptation to environmental changes.'\
\
CCR2_bio = '1. **Extracellular Matrix and Structural Proteins**:\\\
   - **COL1A2, COL3A1, COL6A2, SPARC**: These proteins are involved in the formation and maintenance of the extracellular matrix (ECM). Collagens (COL1A2, COL3A1, COL6A2) provide structural support, while SPARC (Secreted Protein Acidic and Rich in Cysteine) is involved in ECM remodeling and cell-matrix interactions. Their interplay is crucial for tissue integrity and repair processes.\\\
2. **Ribosomal Proteins and Translation**:\\\
   - **RPL23, RPL13, RPLP1, RPS26, TMA7**: These ribosomal proteins are components of the ribosome, essential for protein synthesis. Their coordinated function is critical for maintaining cellular protein homeostasis and responding to increased protein synthesis demands.\\\
3. **Immune Response and Inflammation**:\\\
   - **A2M, IFITM3, CD59, C1S, TAPBP**: A2M (Alpha-2-Macroglobulin) acts as a protease inhibitor and is involved in immune regulation. IFITM3 is known for its role in antiviral defense. CD59 protects cells from complement-mediated lysis, while C1S is part of the complement system. TAPBP (Tapasin) is involved in antigen processing. Together, these proteins suggest a role in modulating immune responses and protecting against pathogens.\\\
4. **Metabolic and Redox Regulation**:\\\
   - **TXNL1, TXNIP, SIRT6, SARDH**: TXNL1 and TXNIP are involved in redox regulation and cellular stress responses. SIRT6 is a sirtuin involved in DNA repair and metabolism, while SARDH is involved in sarcosine metabolism. These proteins may work together to maintain cellular redox balance and metabolic homeostasis.\\\
5. **Cytoskeletal and Cell Adhesion Proteins**:\\\
   - **NES, TUBB8, CFL1, SHANK3, MYZAP**: NES (Nestin) and TUBB8 are involved in cytoskeletal dynamics. CFL1 (Cofilin 1) regulates actin filament dynamics, while SHANK3 and MYZAP are involved in cell adhesion and signaling. These proteins are crucial for maintaining cell structure and facilitating cell-cell communication.\\\
6. **Ubiquitination and Protein Degradation**:\\\
   - **USP17L3, USP17L11, USP17L12, UBR3**: These proteins are involved in ubiquitin-mediated protein degradation. The USP17 family members are deubiquitinating enzymes, while UBR3 is an E3 ubiquitin ligase. They likely regulate protein turnover and cellular protein quality control.\\\
The identified proteins suggest involvement in several key biological processes, including extracellular matrix maintenance, protein synthesis, immune response, redox regulation, cytoskeletal dynamics, and protein degradation. The presence of ECM proteins (COL1A2, COL3A1, COL6A2, SPARC) alongside immune-related proteins (A2M, IFITM3, CD59) may indicate a role in tissue repair and immune modulation. The ribosomal proteins highlight a focus on protein synthesis, potentially in response to cellular stress or increased demand. The interplay between redox regulators (TXNL1, TXNIP) and metabolic enzymes (SIRT6, SARDH) suggests a coordinated response to maintain cellular homeostasis. While no single disease process is immediately apparent, the combination of immune, structural, and metabolic proteins could be relevant in contexts such as inflammation, fibrosis, or cancer, where these processes are often dysregulated.'\
\
Our_mac_process = 'Extracellular matrix organization and immune response modulation\\\
1. **Extracellular Matrix (ECM) Components and Remodeling**: The proteins COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, COL6A1, COL6A3, and BGN are integral components of the ECM, contributing to its structural integrity and function. SPARC and TIMP1 are involved in ECM remodeling, with SPARC modulating cell-matrix interactions and TIMP1 inhibiting matrix metalloproteinases, thus regulating ECM degradation.\\\
2. **Immune Response Modulation**: CIITA is a master regulator of MHC class II gene expression, crucial for antigen presentation and adaptive immune responses. CD74, associated with MHC class II, plays a role in antigen processing. FCGR2C is involved in immune complex binding and clearance. IL10RA is part of the IL-10 receptor complex, mediating anti-inflammatory responses.\\\
3. **Cellular Stress and Apoptosis**: JUNB and JUND are components of the AP-1 transcription factor complex, which regulates gene expression in response to cellular stress and is involved in apoptosis. CFLAR is an inhibitor of apoptosis, modulating cell death pathways.\\\
4. **Protein Synthesis and Ribosomal Function**: RPL19 and RPL10L are ribosomal proteins, essential for protein synthesis. Their presence suggests a role in maintaining cellular protein homeostasis.\\\
5. **Lipid Metabolism and Transport**: ELOVL5 is involved in the elongation of long-chain fatty acids, indicating a role in lipid metabolism.\\\
6. **Signal Transduction and Regulation**: GNB2 and GNAI2 are involved in G-protein signaling pathways, which are critical for transmitting extracellular signals to intracellular responses. DUSP1 is a phosphatase that deactivates MAP kinases, thus regulating signal transduction pathways.\\\
The system of interacting proteins primarily focuses on the structural organization of the ECM and modulation of immune responses, with additional roles in cellular stress response, protein synthesis, and signal transduction. The integration of these processes suggests a coordinated effort to maintain tissue integrity and modulate immune functions.'\
\
Hanxi_mac_bio = 'Collagen and Extracellular Matrix Proteins\\\
- **COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, COL6A1, COL6A3, SPARC, BGN, FN1, SPON2, AEBP1, TIMP1**: These proteins are primarily involved in the formation and maintenance of the extracellular matrix (ECM). Collagens (COL1A1, COL1A2, COL3A1, COL4A1, COL4A2, COL6A1, COL6A3) are structural proteins that provide tensile strength and structural integrity to tissues. SPARC (osteonectin) is involved in collagen binding and ECM remodeling. BGN (biglycan) and FN1 (fibronectin) are ECM components that interact with collagens and other matrix proteins to influence cell adhesion and migration. TIMP1 is a tissue inhibitor of metalloproteinases, regulating ECM degradation. This group suggests a role in tissue remodeling, wound healing, or fibrosis.\\\
Immune System and Inflammation\\\
- **CIITA, IL10RA, CD74, FCGR2C, B2M, CTSZ, LGMN, GRN, CFLAR**: These proteins are associated with immune response and inflammation. CIITA is a transcriptional activator of MHC class II genes, crucial for antigen presentation. IL10RA is a receptor for interleukin-10, an anti-inflammatory cytokine. CD74 is involved in antigen processing and presentation. FCGR2C is a receptor for the Fc region of IgG, playing a role in immune complex clearance. B2M is part of the MHC class I molecule. CTSZ and LGMN are proteases involved in antigen processing. GRN (progranulin) has roles in inflammation and tissue repair. CFLAR is involved in apoptosis regulation. This group suggests involvement in immune regulation and inflammatory processes.\\\
Signal Transduction and Transcription Regulation\\\
- **JUNB, JUND, ZFP36, MAZ, CEBPD, DUSP1, SIRT6, ZNF680, ZNF346, THAP12**: These proteins are involved in signal transduction and transcriptional regulation. JUNB and JUND are components of the AP-1 transcription factor complex, regulating gene expression in response to various stimuli. ZFP36 is involved in mRNA decay, influencing gene expression post-transcriptionally. MAZ is a transcription factor involved in gene regulation. CEBPD is a transcription factor involved in immune and inflammatory responses. DUSP1 is a phosphatase that deactivates MAP kinases, modulating signal transduction. SIRT6 is involved in chromatin remodeling and gene expression regulation. ZNF680, ZNF346, and THAP12 are zinc finger proteins, likely involved in DNA binding and transcriptional regulation. This group suggests roles in cellular response to stress and regulation of gene expression.\\\
Metabolism and Cellular Homeostasis\\\
- **SCO2, ATP5MC2, SAT1, ELOVL5, GSTP1, FTH1, OAZ1**: These proteins are involved in metabolic processes and cellular homeostasis. SCO2 is involved in mitochondrial function and energy production. ATP5MC2 is a component of ATP synthase, crucial for ATP production. SAT1 is involved in polyamine metabolism. ELOVL5 is involved in fatty acid elongation. GSTP1 is a detoxification enzyme. FTH1 is involved in iron storage and homeostasis. OAZ1 regulates polyamine synthesis. This group suggests roles in energy metabolism and cellular stress responses.\\\
The identified proteins suggest a strong involvement in extracellular matrix organization and remodeling, immune response and inflammation, signal transduction and transcriptional regulation, and metabolic processes. The presence of multiple collagens and ECM-related proteins indicates potential roles in tissue remodeling or fibrosis. The immune-related proteins suggest involvement in antigen processing and inflammatory regulation. The transcription factors and signaling molecules indicate regulation of gene expression in response to cellular stress. Metabolic proteins suggest roles in maintaining cellular energy and homeostasis. Collectively, these proteins could be implicated in processes such as wound healing, fibrosis, immune response, and metabolic regulation. There is potential relevance to diseases involving ECM dysregulation, such as fibrosis, and immune-related conditions.'\
\
Our_mac_dist_process = 'Cardiac Muscle Structure and Function\\\
1. **Structural Proteins**: TTN (titin), ACTA1 (actin, alpha skeletal muscle), ACTC1 (actin, alpha cardiac muscle 1), DES (desmin), and MYH7 (myosin heavy chain 7) are integral to the sarcomere, the fundamental unit of muscle contraction. TTN provides elasticity and structural integrity, while ACTA1 and ACTC1 are involved in the contractile apparatus. DES maintains the structural integrity of muscle fibers, and MYH7 is a key component of the thick filament in cardiac muscle.\\\
2. **Regulatory Proteins**: TCAP (titin-cap) interacts with titin and is involved in sarcomere assembly and stability. ANKRD1 (ankyrin repeat domain 1) is a stress response protein that modulates cardiac gene expression in response to mechanical stress.\\\
3. **Signaling and Metabolic Proteins**: AKT3 (protein kinase B gamma) and RPS6KB1 (ribosomal protein S6 kinase beta-1) are involved in signaling pathways that regulate cardiac growth and metabolism. NDUFAB1 (NADH:ubiquinone oxidoreductase subunit AB1) is part of the mitochondrial respiratory chain, crucial for energy production in cardiac cells.\\\
4. **Chaperones and Stress Response Proteins**: HSPA1A (heat shock protein family A member 1A) assists in protein folding and protection against stress-induced damage. TP53INP2 (tumor protein p53 inducible nuclear protein 2) is involved in autophagy, a process important for cellular homeostasis and stress response.\\\
5. **Transcription and Translation Regulation**: MAZ (Myc-associated zinc finger protein) and ELF2 (E74-like factor 2) are transcription factors that regulate gene expression in cardiac cells. EEF1A2 (eukaryotic translation elongation factor 1 alpha 2) is involved in protein synthesis, essential for maintaining cardiac muscle function.\\\
6. **Calcium Handling and Ion Transport**: SLC41A3 (solute carrier family 41 member 3) is involved in magnesium transport, which is critical for cardiac muscle contraction and relaxation.\\\
The proteins in this system are predominantly involved in maintaining the structure, function, and regulation of cardiac muscle, with a focus on sarcomere integrity, signaling pathways, and stress response mechanisms. This integrated network ensures proper cardiac muscle contraction and adaptation to physiological demands.'\
\
Our_mac_dist_bio = '1. **Cardiac and Muscle Function Proteins**: \\\
   - **NPPB, NPPA, TTN, TCAP, ANKRD1, ACTA1, ACTC1, DES, MYH7**: These proteins are primarily associated with cardiac and skeletal muscle function. NPPB and NPPA are natriuretic peptides involved in cardiovascular homeostasis. TTN (titin) and TCAP (telethonin) are structural proteins critical for muscle elasticity and integrity. ANKRD1, ACTA1, ACTC1, DES, and MYH7 are involved in muscle contraction and structural stability. Mutations or dysregulation in these proteins are often linked to cardiomyopathies and muscular dystrophies.\\\
2. **Signal Transduction and Kinase Activity**:\\\
   - **AKT3, MAP2K7, RPS6KB1**: These proteins are involved in signal transduction pathways. AKT3 is part of the PI3K/AKT pathway, which is crucial for cell survival and growth. MAP2K7 is involved in the JNK signaling pathway, which responds to stress stimuli. RPS6KB1 is a kinase that regulates protein synthesis and cell growth. Dysregulation in these pathways is often associated with cancer and metabolic disorders.\\\
3. **Protein Synthesis and Degradation**:\\\
   - **EEF1A2, OTUD4, SNW1, CNOT11**: EEF1A2 is involved in the elongation phase of protein synthesis. OTUD4 is a deubiquitinating enzyme, playing a role in protein degradation. SNW1 is a splicing factor, and CNOT11 is part of the CCR4-NOT complex involved in mRNA degradation. These proteins are crucial for maintaining protein homeostasis, and their dysfunction can lead to neurodegenerative diseases and cancer.\\\
4. **Immune Response and Inflammation**:\\\
   - **IL13, IL32, CD74**: IL13 and IL32 are cytokines involved in immune response and inflammation. CD74 is a receptor that plays a role in antigen presentation. These proteins are often implicated in autoimmune diseases and inflammatory conditions.\\\
5. **Metabolic and Mitochondrial Function**:\\\
   - **NDUFAB1, LRPPRC**: NDUFAB1 is part of the mitochondrial respiratory chain, and LRPPRC is involved in mitochondrial gene expression. Dysfunction in these proteins can lead to metabolic disorders and mitochondrial diseases.\\\
6. **Transcription and Gene Regulation**:\\\
   - **KDM5C, MAZ, ELF2**: KDM5C is a histone demethylase involved in chromatin remodeling. MAZ is a transcription factor, and ELF2 is involved in transcriptional regulation. Alterations in these proteins can affect gene expression and are linked to intellectual disabilities and cancer.\\\
The proteins identified are primarily involved in cardiac and muscle function, signal transduction, protein synthesis and degradation, immune response, metabolic processes, and transcriptional regulation. The strong presence of muscle-related proteins suggests a potential focus on muscle integrity and function, possibly indicating a role in cardiomyopathies or muscular dystrophies. The involvement of signaling and transcriptional regulation proteins points to potential links with cancer and metabolic disorders. Immune-related proteins suggest a possible connection to inflammatory diseases. Overall, the interplay of these proteins suggests a complex network that could be involved in maintaining cellular homeostasis, with implications for various disease processes, particularly those affecting the heart, muscles, and immune system.'\
\
\
comp_wordcloud_bio = create_comparative_wordcloud_color(\
    preprocess_terms(L_ICM_fib_process).join(preprocess_terms(L_ICM_fib_bio)),\
    preprocess_terms(AF_process).join(preprocess_terms(AF_bio)),\
    'Lavine ICM v CTRL Fibroblasts', 'AF Fibroblasts'\
    )\
# Save the image\
plt.savefig("comparative_wordcloud_L_ICM_v_AF_fib.png", dpi=300, bbox_inches='tight')}