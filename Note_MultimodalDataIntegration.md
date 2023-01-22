## Summary
To search for patients with non-small-cell lung cancer who are most likely to respond to immunotherapy, authors develop a computational workflow to extract patient-level features from the multimodal baseline data including **expert-annotated** **computed tomography (CT) scan images**, digitized programmed death-1 (PD-1) and programmed death ligand-1 (PD-L1) **immunohistochemistry (IHC) slides** and **genomic sequencing data**, and outcomes of immune checkpoint inhibitor therapy, and then use a deep learning technique to integrate multimodal features into a **multimodal risk prediction model**.
  - **Data**  
    - To characterize the spatial distributions of CT signal intensities for each tumor, authors extract the quantitative image features from the segmented areas.
    - To characterize the spatial organization of PD-L1 immunoreactivity within the tumor, authors use PD-L1 immunohistochemistry slides.
    - The modes of data are one-dimensional feature vectors derived from PD-L1 IHC, segmented CT scans, and genomic alterations. 
    - **Limitation**: 1. insufficient sample size limit the generalizability of results. 2. small sample size also leads to data highly processed using expert knowledge, which hinders the model to learn the best features associated with immunotherapy reponse in a data-driven manner.
  - **Model**
    - The author's deep learning model stems from an article of **attention-based deep multiple instance learning**. 
    - The strength of this model: 1. it does not restrict the number of data types, enabling the model to be extended to other data modalities and to be applied to other cancers and diseases. 2. it considers the modalities for each patient with **different weights**, maximizating the predictive capacity of each type of data and thus improving the accuracy in predicting the response to immunotherapy (vs single biomarker). 3. even a patient has missing values in the CT scan or genomic sequence, the model can mask the missing data modality and compute a risk prediction. 4. the model can classify patients as immunotherapy responders or non-responders early after treatment.
    - **Performance**: as the table shown below 
    
    
        | Model                                       | AUC           | 95% CI        |
        | -------------                               | ------------- | ------------- |
        | Multimodal                                  | 0.80          | 0.74-0.86     |
        | Tumor mutational burden (unimodal)          | 0.61          | 0.52-0.70     |
        | PD-L1 immunohistochemistry score (unimodal) | 0.73          | 0.65-0.81     |


## Technical details
- **Area under the curve (AUC)**: This quantity is used to measure the performance of a biomarker on distinguishing response versus nonresponse. In general, AUC or Receiver operating characteristic (ROC) curve compares and evaluates the performance of a binary classification model. The curve can be expressed in a plot of **sensitivity/true positive rate** (Y-axis) versus **1-specificity/false positive rate** (X-axis) at different probability cutoffs. The higher AUC suggests the better performance of the classifier. The diagonal line represents the random classification model, where all points along the diagonal line suggest the same true positive and false positive rate.

    | Cut-off	                | Sensitivity               | Specificity               | 1 - Specificity           |
    |-------------------------|---------------------------|---------------------------|---------------------------|
    | 0                       | 1                         | 0                         | 1                         |
    | 0.01                    | 0.979                     | 0.081                     | 0.919                     |
    | 0.02                    | 0.938                     | 0.158                     | 0.842                     |  
    | ...                     | ...                       | ...                       | ...                       |
    | 0.99                    | 0.02                      | 0.996                     | 0.004                     |
    | 1                       | 0                         | 1                         | 0                         |

- **Permutation-tested AUC**: 
  - Context: assume two classifiers using the same test dataset are created and each one has its own distinct ROC curve and AUC value. The classifier M shows a higher AUC value than the classifer N. 
  - Question: is this difference in the AUC values systematic, or random? 
  - Statistical translation: what is the probability of observing this difference or more extreme difference under the null hypothesis? 
  - Computational algorithm using the randomized permuation: 
    - for i = 0, ..., n - 1
      - for cutoff = 0, ..., 1
        - M = (1 - specificity, sensitivity)_cutoff + M
        - N = (1 - specificity, sensitivity)_cutoff + N
      - Z_i = a random pairwise mix between M and N
      - D_i = AUC of Z_i - AUC of N
      - if D_i >= AUC of M - AUC of N
        - then count = count + 1
     - p = count/n
     - If p <= threshold (e.g., 0.05), we may conclude the difference in the AUC values of two classifiers is statistically significant.
  
- Repeated subsampling-tested AUC: AUC calculated using the test data from the repeated subsampling. Random subsampling, also known as Monte Carlo crossvalidation, refers to randomly splitting the data into subsets repeatedly, where the users define the size of the subsets. In contrast to a traditional cross-validation procedure, its advantage is resulting in more pessimistic predictions of the test data compared with cross-validation. The predictions give a realistic estimation of the predictions of external validation data.
- Autocorrelation feature
- **Gray level co-occurrence matrix (GLCM) features**: GLCM is defined over an image to be the distribution of co-occurring pixel values (grayscale values, or colors) at a given offset. Authors exploit GLCMs, commonly used in image processing to quantify the similarity of neighboring pixels, to characterize PD-L1 expression.
- **Multiple instance logistic regression (MILR)**：MILR is a machine learning approach where sets of training data (bags) share a common label. Authors decorates MILR with the attention-based pooling, a strategy that assigns attention-based weights to each instance. These weights are a function of the instance features and are optimized for the outcome prediction. They are dynamically determined. The attention-based pooling is applicable for the fixed number of features while the unfixed number of input instances. The final output score is a weighted sum of the same logistic regression model for each instance. MILR treats each instance equally as the single set of parameters are shared across all instances. This machine learning technique has some **limitations**. First of all, it typically handles a single bag of homogenous instances, so heterogenous instances are not suitable for this technique. Second, the single set of parameters shared across all instances may not capture the possible difference in the relationship between the lesions of the same shape in different sites and the predictor.
  - Model parameters: hidden size: 32; balanced class weights; binary cross-entropy loss; 0.005 learning rate; 250 steps; 0.005 L2 regularization strength using the Adam optimizer.
- **Dynamic attention with masking (DyAM)**: Authors assign attention-based weights to each input data mode so that each mode has its own 2 sets of trainable parameters: 1. risk parameters that map the input vectors to the label; 2. attention parameters that map the input vector to an attention score. The final output score is a weighted sum of LR models applied for each input mode. The **advantages** of DyAM, compared to MILR, are: 1. because attention-based weights are trainable and DyAm treats each mode specifically while still optimizing all parameters in concert, the model learns which input modes are most relevant to treatment response; 2. To address missing data (e.g., non-segmentable disease in CT), DyAM adds a masking function that sets any missing mode's attention weights to zero.
  - Model parameters: balanced class weights; binary cross-entropy loss; 0.01 learning rate; 125 steps; 0.001 L2 regularization strength using the Adam optimizer.
- **Tumor proportion score (TPS)**: this measurement is defined as the percent of partial or complete membranous staining among viable tumor cells. A negative score refers to the staining in <1% of tumor cells or the absence of staining in tumor cells. Authors exclude the PD-L1 IHC slides that did not meet the minimum number of tumor cells for PD-L1 TPS assessment (<100 tumor cells).
- **DenseNet AI V2 classifier**：DenseNet is a deep learning technique to classify images. It has a network architecture where each layer is directly connected to every other layer in a feed-forward fashion (within each dense block). For each layer, the feature maps of all preceding layers are treated as separate inputs whereas its own feature maps are passed on as inputs to all subsequent layers.
