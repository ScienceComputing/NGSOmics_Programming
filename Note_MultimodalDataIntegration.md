## Summary
To search for patients with non-small-cell lung cancer who are most likely to respond to immunotherapy, authors develop a computational workflow to extract patient-level features from the multimodal baseline data including **expert-annotated** **computed tomography (CT) scan images**, digitized programmed death-1 (PD-1) and programmed death ligand-1 (PD-L1) **immunohistochemistry (IHC) slides** and **genomic sequencing data**, and outcomes of immune checkpoint inhibitor therapy, and then use a deep learning technique to integrate multimodal features into a **multimodal risk prediction model**.
  - **Data**  
    - To characterize the spatial distributions of CT signal intensities for each tumor, authors extract the quantitative image features from the segmented areas.
    - To characterize the spatial organization of PD-L1 immunoreactivity within the tumor, authors use PD-L1 immunohistochemistry slides.
    - **Limitation**: 1. insufficient sample size limit the generalizability of results. 2. small sample size also leads to data highly processed using expert knowledge, which hinders the model to learn the best features associated with immunotherapy reponse in a data-driven manner.
  - **Model**
    - The author's deep learning model stems from an article of **attention-based deep multiple instance learning**. 
    - The strength of this model: 1. it does not restrict the number of data types, enabling the model to be extended to other data modalities and to be applied to other cancers and diseases. 2. it considers the modalities for each patient with **different weights**, maximizating the predictive capacity of each type of data and thus improving the accuracy in predicting the response to immunotherapy (vs single biomarker). 3. even a patient has missing values in the CT scan or genomic sequence, the model can mask the missing data modality and compute a risk prediction. 4. the model can classify patients as immunotherapy responders or non-responders early after treatment.
    - **Performance**: as shown below 
    
    
    | Model                                       | AUC           | 95% CI        |
    | -------------                               | ------------- | ------------- |
    | Multimodal                                  | 0.80          | 0.74-0.86     |
    | Tumor mutational burden (unimodal)          | 0.61          | 0.52-0.70     |
    | PD-L1 immunohistochemistry score (unimodal) | 0.73          | 0.65-0.81     |


## Terminology
- **Dynamic attention with masking (DyAM)**
- Permutation-tested area under the curve (AUC)
- Repeated subsampling-tested AUC
- Autocorrelation feature
- **Gray level co-occurrence matrix (GLCM) features**: GLCM is defined over an image to be the distribution of co-occurring pixel values (grayscale values, or colors) at a given offset. Authors exploit GLCMs, commonly used in image processing to quantify the similarity of neighboring pixels, to characterize PD-L1 expression.
- **Multiple instance logistic regression (MILR)**：MILR is a machine learning approach where sets of training data (bags) share a common label. Authors decorates MILR with the attention-based pooling, which assign attention-based weights to each instance. These weights are a function of the instance features and are optimized for the outcome prediction. They are dynamically determined. The attention-based pooling is applicable for the fixed number of features while the unfixed number of input instances. The final output score is a weighted sum of the samle logistic regression model for each instance. MILR treats each instance equallly as the single set of parameters are shared across all instances.
- **Tumor proportion score (TPS)**: this measurement is defined as the percent of partial or complete membranous staining among viable tumor cells. A negative score refers to the staining in <1% of tumor cells or the absence of staining in tumor cells. Authors exclude the PD-L1 IHC slides that did not meet the minimum number of tumor cells for PD-L1 TPS assessment (<100 tumor cells).
- **DenseNet AI V2 classifier**：DenseNet is a artificial intelligence technique to classify images. It has a network architecture where each layer is directly connected to every other layer in a feed-forward fashion (within each dense block). For each layer, the feature maps of all preceding layers are treated as separate inputs whereas its own feature maps are passed on as inputs to all subsequent layers.
