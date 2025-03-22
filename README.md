# Gene Expression Based Subtyping and Prognostic Analysis in Breast Cancer
This study focuses on classifying breast cancer subtypes using gene expression data and identifying key differentially expressed genes. Biological significance is examined to understand their role in cancer progression. The analysis also explores how breast cancer subtypes influence survival outcomes in relation to clinical attributes.

# Objectives

1. Classify Breast Cancer Subtypes based on gene expression profiles.

2. Identify Key Differentially Expressed Genes that are statistically and biologically significant across subtypes.

3. Investigate the Biological Role of key genes in breast cancer progression.

4. Assess Prognostic Impact, evaluating how breast cancer subtypes influence survival outcomes in relation to clinical attributes.


# Background
Breast cancer encompasses a group of diseases traditionally classified using clinical and pathological criteria, such as histological grade and lymph node metastasis status (LN status). Additionally, the expression of receptor proteins, including estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor receptor 2 (HER2), informs treatment strategies and prognosis (Clarke et al., 2013).

Advancements in translational ‘omics’ have addressed breast cancer heterogeneity through genomic, transcriptomic, and proteomic profiling, leading to the identification of distinct molecular subgroups (Curtis et al., 2012; Kao et al., 2011). These subgroups exhibit unique biological characteristics and are associated with varying clinical outcomes (Asleh et al., 2022; Van ’T Veer et al., 2002). Refined classification based on gene expression offers valuable insights and holds potential for improving personalised breast cancer treatment.


# Methods  

## Data Description  
The dataset consists of **single-channel microarray gene expression data** covering **45,685 genomic locations** from **251 breast cancer patients**. Clinical features and survival data were also provided.  

All analyses and visualisations were performed using **Bioconductor v3.17** (**BiocManager v1.30.22**) and **R v4.3.1 (2023-06-16)**.  

## Data Pre-processing and Expression Cluster Analysis  

- **Expression Matrix Construction**  
  - Imported data files and inspected them.  
  - Created an **expression matrix** with gene names as row indices using **Affymetrix Human Genome U133 Plus 2.0 Array** annotations.  
  - Performed **between-array normalisation** to ensure comparability of gene expression intensities across samples.  
  - Scaled the expression data to standardise gene-level variance.  

- **Clustering Analysis**  
  - Computed a **distance matrix** using the **Euclidean method**.  
  - Performed **hierarchical clustering** with the **Ward.D method**.  
  - Used the **silhouette function** to determine the optimal number of clusters.  
  - Conducted **Principal Component Analysis (PCA)** using the **top 10% (4,500 genes)** with the highest variance.  
  - Generated a **heatmap** using the **top 1,500 most variable genes** to explore cluster relationships.  

## Differential Gene Expression (DE) Analysis  

- Identified **differentially expressed genes (DEGs)** between clusters using the **limma** package.  
- Applied the **Empirical Bayes method** to the DE object.  
- Generated a **volcano plot** to visualise differentially expressed genes.  
- Applied **False Discovery Rate (FDR) correction** at **0.05**.  
- Retrieved genes with **q-values ≤ 0.05**.  

## Enrichment Analysis  

- Extracted the **top 100 differentially expressed genes** for **Gene Ontology (GO) enrichment analysis**.  
- Conducted **GO analysis** using the **Gene Ontology Resource website** for biological processes in **Homo sapiens**.  
  - Selected **“PANTHER GO-Slim Biological Process”** as the annotation dataset.  
  - No correction was applied.  
- Performed **GO annotation analysis** using **GOnet** (Pomaznoy et al., 2018) to generate interactive gene-GO term graphs.  
  - **Species:** Human  
  - **GO Namespace:** “biological_process”  
  - **Analysis Type:** “GO term annotation”  
  - **Annotation Analysis:** “Generic GO slim”  

## Survival Analysis  

- Computed **gene scores** by summing and normalising expression values of **significant DEGs**.  
- Used a **Cox regression model** to calculate **hazard ratios** between expression clusters and survival data.  
- Incorporated **clinical features** into the model to assess their impact.  
- Generated **Kaplan-Meier survival curves** to visualise different survival outcomes across clusters.  

# Results

## Cluster Identification
- The optimal number of clusters, maximizing the median silhouette value, was determined to be 3, based on analysis of cluster numbers between 2 and 5 (Figure 1A).
- Hierarchical clustering (Ward’s method) identified three distinct subtypes.
- Principal Component Analysis (PCA) of genes with the highest expression variance, as shown in Figure 1B, illustrates the similarity between the identified subtypes in 2-dimensional space.

**A** 

![fig1A](https://github.com/user-attachments/assets/30071cb8-de81-47b2-98a3-4dd57ab3967e)

**B**

![fig1B](https://github.com/user-attachments/assets/a2a63a00-378f-44ca-8a55-4ab3be0fd72a)


Figure 1. Cluster analysis and diagnostics. (A) silhouette plot for breast cancer patients’ expression data. (B) Principal component analysis for expression data.


## Heatmap
- A heatmap generated from the scaled expression data, using the top 1500 genes with the highest variance, shows cluster memberships determined by Ward’s method (Figure 2).
- Cluster C, represented as downregulated (in red hues), is underexpressed in cluster 1 compared to clusters 2 and 3.

![Screenshot 2025-03-20 at 21 08 44](https://github.com/user-attachments/assets/6e269c0f-62ca-4404-a45f-9770cbb10b81)



Figure 2. Relationship between clusters generated with high variance expression data and patient clusters identified by Ward's method.

## Differential Gene Expression (DE)
- The volcano plot generated by the **limma** package shows that most differentially expressed genes (DEGs) are highly statistically significant (log(p-value) > 50) and upregulated (logFC > 2) in clusters 2 and 3 when compared with cluster 1 (Figure 3).
- After applying FDR correction at a 0.05 significance level:
  - 22,435 genes were differentially expressed between clusters 1 and 2 (p-value ≤ 0.05) and 23,714 genes with q-value ≤ 0.05.
  - Between clusters 1 and 3, 22,269 genes were differentially expressed (p-value ≤ 0.05) and 24,084 genes with q-value ≤ 0.05.
- Notably, more than 50% of the genomic locations are overexpressed across the studied samples.

![image](https://github.com/user-attachments/assets/2ddbe5ed-5948-4503-80e1-ee9b8235907b)

Figure 3. Volcano plot of differentially expressed genes between cluster 1 and cluster 2 and 3 combined.



## Enrichment Analysis
- The **PANTHER** gene ontology (GO) enrichment analysis of the top 100 DEGs, shown in Figure 4, highlights significant involvement in **cell cycle processes**, including **DNA damage checkpoint signalling** and **negative regulation of cell cycle**.
- A notable proportion of patients exhibited **overexpression of estrogen (ER)** and **progesterone (PR) receptors**, with **84.86%** and **75.7%** of patients, respectively, showing overexpression.


![image](https://github.com/user-attachments/assets/5dad36e5-7d48-489f-930a-30d20d2097d3)


Figure 4.  PANTHER multiple pie chart view generated by Gene ontology resource website for top 100 DE genes with GO-Slim Biological Process




## GO Terms and Network Analysis
- Figure 5 illustrates the key functions related to **cell growth**, **cell cycle**, and **proliferation**, with a focus on **ESR1**, the estrogen receptor gene, which is predominantly overexpressed in ER-positive patients.

![image](https://github.com/user-attachments/assets/055d5688-158e-4c82-a43a-8868cac95f4c)

Figure 5. Exploration of genes and GO terms as a graph generated by GOnet (a tool for interactive Gene Ontology analysis).



## Survival Analysis
- The **Kaplan-Meier survival curve** (Figure 6) shows that survival probability is lower for patients in cluster 1 compared to clusters 2 and 3.
- The **Cox regression model** calculates hazard ratios of 1.4284 for cluster 1 (>1), indicating increased risk, while clusters 2 and 3 have hazard ratios of 0.3947 and 0.5172 (<1), respectively, suggesting a protective effect from their gene expression.

![fig6](https://github.com/user-attachments/assets/93235ce3-3103-4a5d-8c85-0859418eb130) ![image](https://github.com/user-attachments/assets/5d14fe2c-d953-420f-91b6-974a1df8ad86)

Figure 6. Survival curve for breast cancer patients stratified by cluster membership (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.




- **Lymph node metastasis status (LN status)** also significantly affects survival, with LN-positive patients exhibiting lower survival rates (Figure 7).

![fig7](https://github.com/user-attachments/assets/47d5a353-11c6-463e-a300-560f19a3f4cd)![image](https://github.com/user-attachments/assets/9730eb1d-3aa6-4513-9833-a655e1249a3f)

Figure 7. Survival curve for breast cancer patients stratified by LN status (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.



- When combined with cluster membership, LN-positive patients in cluster 1 have significantly higher hazard ratios (Figure 8).


![fig8](https://github.com/user-attachments/assets/cd77a10b-bc01-4335-abe6-1a71f8b4846b) 
![image](https://github.com/user-attachments/assets/0f84182d-8fed-4dc8-9e6d-b55777e144ba)

Figure 8. Survival curve for breast cancer patients stratified by LN status and cluster membership (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.



### Tumor Size
- Survival analysis by tumor size (≥ 20 mm vs < 20 mm) shows that patients with larger tumors have a poorer prognosis, with significantly higher hazard ratios for this group (Figure 9).
![fig9](https://github.com/user-attachments/assets/75fecadc-cfbd-47bc-9d14-2760b1e4db05)

![image](https://github.com/user-attachments/assets/598c537e-7597-432c-9d40-5d45717ab2ab)

Figure 9. Survival curve for breast cancer patients stratified by Tumour size > 20 mm (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.

- Moreover, patients in cluster 1 and those with larger tumors exhibit the lowest survival rates (Figure 10).

![fig10](https://github.com/user-attachments/assets/b724e9f3-18a7-4a3d-9729-097558f42085)

![image](https://github.com/user-attachments/assets/bf48fc61-fbf1-4f38-8e21-642553e3dab4)

Figure 10. Survival curve for breast cancer patients stratified by Tumour size > 20 mm and cluster membership (top) Cox regression summary table with hazard ratios (exp(coef)) (bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.'




### Age at Diagnosis
- The survival analysis based on age at diagnosis did not show significant results, with the median age of diagnosis being around **62 years**.
- There was no evidence suggesting an improved prognosis for those diagnosed at an earlier age (Figure 11).


![image](https://github.com/user-attachments/assets/f0ef9b1f-cbd4-43ff-a889-75b58ad5d28f)

![image](https://github.com/user-attachments/assets/b138bce0-b35c-41d7-af5b-8c2227ea513f)

Figure 11. Survival curve for breast cancer patients stratified by Age at diagnosis (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.


### ER and PR Status
- Survival analysis for **ER** and **PR** status shows no significant differences, possibly due to imbalanced data. However, as previously mentioned, a higher proportion of patients showed overexpression of these receptor proteins (Figures 12 & 13).

![image](https://github.com/user-attachments/assets/2c842c61-4f2f-492a-85c3-68b1d0853be3)
![image](https://github.com/user-attachments/assets/9f023c5c-f697-4b7a-8c1d-1f32e54e1f40)

Figure 12. Survival curve for breast cancer patients stratified by ER expression status (top) Cox regression summary table with hazard ratios (exp(coef)(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.


![image](https://github.com/user-attachments/assets/1a861df5-823e-4177-b39e-c3d65d7bd734)
![image](https://github.com/user-attachments/assets/f74053e0-2205-4c5d-8d17-457b752a6d2c)

Figure 13. Survival curve for breast cancer patients stratified by PgR status (top) Cox regression summary table with hazard ratios (exp(coef)) (bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.



### Histological Grade
- Kaplan-Meier survival curves stratified by histological grade (Figure 14) show that **G1** patients have the highest survival rates, while **G3** patients exhibit the lowest.
- However, these results could not be explained by gene expression effects, as the hazard ratios were insignificant.

![image](https://github.com/user-attachments/assets/c0c14e87-b0c1-4edf-97d5-a563380f5077)
![image](https://github.com/user-attachments/assets/c3ede2f3-ad2e-4d7b-805c-aa1ba194ac0c)

Figure 14. Survival curve for breast cancer patients stratified by Histologic grade (top) Cox regression summary table with hazard ratios (exp(coef))(bottom). Significant codes for p-values are indicated as:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1.

# Conclusion

Hierarchical clustering of microarray transcription data identified three subtypes among the breast cancer patients in this dataset. Clusters 2 and 3 are significantly different from cluster 1 in both gene expression and survival outcomes. The differential expression analysis, combined with Gene Ontology (GO) enrichment, suggests that members of clusters 2 and 3 show upregulation of genes associated with cell cycle and proliferation compared to cluster 1. This finding is further supported by the heatmap, which highlights the differential upregulation of genes in clusters 2 and 3. Gene expression profiles of cluster 1 are associated with poorer survival rates compared to the other two clusters, with cluster 2 showing the highest protective effect. When gene expression levels are combined with other clinical variables, such as lymph node metastasis status (LN status) and tumor size greater than 20 mm, it is evident that patients with these characteristics exhibit significantly higher hazard ratios.

# Discussion

Various methods were explored to identify subgroups in the dataset, with all available algorithms in the `hclust` {stats} function being tested sequentially. The silhouette method was used to determine the optimal number of clusters, and principal component analysis (PCA) was performed to visualize the separation of clusters in 2-dimensional space. Based on these analyses, the "ward.D" method was chosen for hierarchical clustering.

Three tools were employed in GO enrichment analysis. In addition to the methods described in the methodology, the **Gene Ontology Analysis** using the `goana` {limma} function was applied to the list of differentially expressed (DE) genes, followed by the **Top GO Terms** function from `topGO` {limma}. This method highlighted several biological processes, including metabolic pathways, but no significant conclusions could be drawn from them. Therefore, the top 100 Entrez Gene IDs, extracted from the `topTable` {limma} function applied to the linear model fit produced by `lmFit` {limma}, were used in both GO analysis processes. For GO analysis, the "GO-Slim Biological Process" was chosen over "GO complete Biological Process" as it offers a concise, manually curated set of biological processes.

The overrepresentation of GO terms related to the cell cycle and proliferation, including negative regulation of the cell cycle, strongly suggests a tumorous trait in the dataset, as aberrant activation of these processes is a hallmark of cancer cells (Hanahan & Weinberg, 2011). Additionally, the upregulation of **ESR1**, implicated in estrogen receptor (ER) signaling pathways, was observed in clusters 2 and 3, as compared to cluster 1. When examining the frequency of **ER+** status, cluster 1 exhibited a lower proportion of ER+ patients than clusters 2 and 3 (see supplementary R code). Moreover, the **Luminal/ER+** molecular subtype identified by Hu et al. (2006) lists **ER**, **FOXA1**, **GATA3**, and **XBP1** as key markers, and three of these four genes (ER, FOXA1, GATA3) were found in the top DE genes of clusters 2 and 3 (Figure 5). **GATA3** was also identified as a key gene in the **Luminal A subtype** by Jia et al. (2020), which is clinically characterized by higher survival rates (Orrantia-Borunda et al., 2022). This is consistent with the survival analysis results of clusters 2 and 3. However, the literature linking these findings is not comprehensive, and further review is needed.

Two other genes, **CDK1** and **IDO1**, were present in the top DE gene list. The overexpression of these genes has been strongly associated with breast cancer progression (Feng et al., 2020; Sofi et al., 2022). Additionally, **EGR1** was overexpressed in the DE gene list. Recent reviews have recognized **EGR1** as both a tumor suppressor and an oncogene in various cancers, with conflicting evidence regarding its role in cancer progression (Wang et al., 2021). Further studies are needed to clarify the mechanism of **EGR1** in cancer initiation and progression.


# References:
Asleh, K., Riaz, N., & Nielsen, T. O. (2022). Heterogeneity of triple negative breast cancer: Current advances in subtyping and treatment implications. Journal of Experimental & Clinical Cancer Research, 41(1), 265. https://doi.org/10.1186/s13046-022-02476-1

Clarke, C., Madden, S. F., Doolan, P., Aherne, S. T., Joyce, H., O’Driscoll, L., Gallagher, W. M., Hennessy, B. T., Moriarty, M., Crown, J., Kennedy, S., & Clynes, M. (2013). Correlating transcriptional networks to breast cancer survival: A large-scale coexpression analysis. Carcinogenesis, 34(10), 2300–2308. https://doi.org/10.1093/carcin/bgt208

Curtis, C., Shah, S. P., Chin, S.-F., Turashvili, G., Rueda, O. M., Dunning, M. J., Speed, D., Lynch, A. G., Samarajiwa, S., Yuan, Y., Gräf, S., Ha, G., Haffari, G., Bashashati, A., Russell, R., McKinney, S., Langerød, A., Green, A., Provenzano, E., ... Aparicio, S. (2012). The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature, 486(7403), 346–352. https://doi.org/10.1038/nature10983
13

Feng, X., Tang, R., Zhang, R., Wang, H., Ji, Z., Shao, Y., Wang, S., Zhong, T., Gu, Y., & Meng, J. (2020). A comprehensive analysis of IDO1 expression with tumour‐ infiltrating immune cells and mutation burden in gynaecologic and breast cancers. Journal of Cellular and Molecular Medicine, 24(9), 5238–5248. https://doi.org/10.1111/jcmm.15176

Gene Ontology Resource. (n.d.). Retrieved October 24, 2023, from https://geneontology.org/ 

GOnet. (n.d.). Retrieved October 23, 2023, from https://tools.dice-database.org/GOnet/ 

Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of Cancer: The Next Generation. Cell,
144(5), 646–674. https://doi.org/10.1016/j.cell.2011.02.013

Hu, Z., Fan, C., Oh, D. S., Marron, J., He, X., Qaqish, B. F., Livasy, C., Carey, L. A.,
Reynolds, E., Dressler, L., Nobel, A., Parker, J., Ewend, M. G., Sawyer, L. R., Wu, J., Liu, Y., Nanda, R., Tretiakova, M., Orrico, A. R., ... Perou, C. M. (2006). The molecular portraits of breast tumors are conserved across microarray platforms. BMC Genomics, 7(1), 96. https://doi.org/10.1186/1471-2164-7-96

Jia, R., Li, Z., Liang, W., Ji, Y., Weng, Y., Liang, Y., & Ning, P. (2020). Identification of key genes unique to the luminal a and basal-like breast cancer subtypes via bioinformatic analysis. World Journal of Surgical Oncology, 18(1), 268. https://doi.org/10.1186/s12957-020-02042-z

Kao, K.-J., Chang, K.-M., Hsu, H.-C., & Huang, A. T. (2011). Correlation of microarray- based breast cancer molecular subtypes and clinical outcomes: Implications for treatment optimization. BMC Cancer, 11(1), 143. https://doi.org/10.1186/1471-2407- 11-143

Orrantia-Borunda, E., Anchondo-Nuñez, P., Acuña-Aguilar, L. E., Gómez-Valles, F. O., & Ramírez-Valdespino, C. A. (2022). Subtypes of Breast Cancer. In H. N. Mayrovitz (Ed.), Breast Cancer. Exon Publications.
http://www.ncbi.nlm.nih.gov/books/NBK583808/

Pomaznoy, M., Ha, B., & Peters, B. (2018). GOnet: A tool for interactive Gene Ontology
analysis. BMC Bioinformatics, 19(1), 470. https://doi.org/10.1186/s12859-018-2533-3 

Sofi, S., Mehraj, U., Qayoom, H., Aisha, S., Almilaibary, A., Alkhanani, M., & Mir, M. A. (2022). Targeting cyclin-dependent kinase 1 (CDK1) in cancer: Molecular docking
and dynamic simulations of potential CDK1 inhibitors. Medical Oncology (Northwood, London, England), 39(9), 133. https://doi.org/10.1007/s12032-022- 01748-2

Van ’T Veer, L. J., Dai, H., Van De Vijver, M. J., He, Y. D., Hart, A. A. M., Mao, M., Peterse, H. L., Van Der Kooy, K., Marton, M. J., Witteveen, A. T., Schreiber, G. J., Kerkhoven, R. M., Roberts, C., Linsley, P. S., Bernards, R., & Friend, S. H. (2002). Gene expression profiling predicts clinical outcome of breast cancer. Nature, 415(6871), 530–536. https://doi.org/10.1038/415530a

Wang, B., Guo, H., Yu, H., Chen, Y., Xu, H., & Zhao, G. (2021). The Role of the Transcription Factor EGR1 in Cancer. Frontiers in Oncology, 11. https://www.frontiersin.org/articles/10.3389/fonc.2021.642547
