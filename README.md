# ConCon_Alignment

## Goal: 
This Matlab code aligns cortical surfaces based on continuous connectivity (ConCon) profiles. 
Example structural connectivity data, stored as matrices, can be downloaded from:

 **[Example Data][https://www.dropbox.com/scl/fo/wpkt84ik2x4wdw3rswqc0/AIkUEPTum07a-ZEV5cmYzfM?rlkey=a56mtlcria97oi26fw5362430&st=uwrax1vi&dl=0]**
 
We provide data for 10 example subjects, each presented as a $4121\times 4121$ upper triangle matrix. Since the structual connectivity is symmetric, you can recover the whole matrix by summing the matrix and its transpose.
## Project Structure  

| Folder        | Description |
|--------------|-------------------------------------------------------------|
| `example/`   | Stores example data (downloaded from Dropbox) |
| `grid/`      | Contains left and right surface grids |
| `results/`   | Stores results, including templates and warping files |
| `scripts/`   | Contains all processing scripts |
| `subject_ids/` | Stores the file with all HCP subject IDs |


## How to Run
### Manually Execution
1. **Download Example Data**  
   - Download the dataset from the provided Dropbox link.  
   - Place the files in the `./example` directory.  
2. **Generate a Template**  
   - Run the following script to create a template from the example data:  
     ```matlab
     run('./scripts/pop_template.m')
     ```
   - The generated `template.mat` file will be saved in the `./results` directory.  
3. **Align Subjects to the Template**  
   - Execute the registration script to align all example subjects to the generated template:  
     ```matlab
     run('./scripts/register_test.m')
     ```
4. **PCA Test**  
   - Run the PCA test script to compute PCA scores for both:  
     1. Concatenated, flattened original structural matrices  
     2. Concatenated, flattened registered structural matrices  
     3. Concatenated, geodesic distance between the original mesh and the warped mesh vertices

   - Execute the script using:  
     ```matlab
     run('./scripts/pca_test.m')
     ```


