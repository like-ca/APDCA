# APDCA
Package: APDCA
Type: Package
Title: APDCA: A novel computational method for predicting associations between Long non-coding RNAs and diseases
Version: 1.0
Author: YangsongHe <yangsonghe03@gmail.com>
Maintainer: YangsongHe <yangsonghe03@gmail.com>
Description: APDCA is an algorithm designed for predicting long non-coding RNA (lncRNA)-disease associations using matrix factorization with enforced sparsity constraints. The method leverages a Difference of Convex (DC) functions approach to incorporate sparsity within a matrix tri-factorization framework.
        
Depends:
    MATLAB (>= 2012a)
License: All source code is copyright, under the Artistic-2.0 License.
		For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)
		
To operate:
	* for the example of APDCA:
		we have provided some pre-data for APDCA, if you want to run APDCA framwork, please run the script "APDCA_Init" directly.      
    		
	**  All data loading processes are included in script "APDCA_Init" 

Files:

APDCA_Init.m: The entry to the APDCA framwork, for loading data and doing some preprocessing.

APDCA_Demo.m: The main function for the algorithm.

- If you have any problem, please contact YangsongHe <yangsonghe03@gmail.com>!
