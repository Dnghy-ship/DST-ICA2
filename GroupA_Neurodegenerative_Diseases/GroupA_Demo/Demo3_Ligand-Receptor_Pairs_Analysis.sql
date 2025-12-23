/*
Overall Function:
Get ligand-receptor pair data for a single disease in a specific tissue
Returns key info (expression levels, significance) of ligand-receptor pairs linked to the disease/tissue
*/

DROP PROCEDURE IF EXISTS GetDiseaseLigandReceptor;

DELIMITER $$

CREATE PROCEDURE GetDiseaseLigandReceptor(
    IN disease_code VARCHAR(10),
    IN tissue_name VARCHAR(100),
    IN max_results INT
)
BEGIN
    DECLARE condition_id_val INT;
    
    
    SELECT condition_id INTO condition_id_val 
    FROM Condition_table 
    WHERE disease = disease_code AND tissue = tissue_name;
    
    IF condition_id_val IS NULL THEN
        SELECT CONCAT('No data found for ', disease_code, ' in ', tissue_name) AS message;
    ELSE
        
        SELECT 
            lr.sig_pair_id,
            lr.lr_pair,
            lr.pattern,
            lg.symbol AS ligand_symbol,
            lg.description AS ligand_description,
            rg.symbol AS receptor_symbol,
            rg.description AS receptor_description,
            ROUND(lr.ligand_log2FC, 3) AS ligand_log2FC,
            ROUND(lr.ligand_pvalue, 6) AS ligand_pvalue,
            ROUND(lr.receptor_log2FC, 3) AS receptor_log2FC,
            ROUND(lr.receptor_pvalue, 6) AS receptor_pvalue
        FROM Ligand_receptor_table lr
        JOIN Gene_table lg ON lr.ligand_GNO = lg.GNO
        JOIN Gene_table rg ON lr.receptor_id = rg.GNO
        WHERE lr.condition_GNO = condition_id_val
        ORDER BY (lr.ligand_pvalue + lr.receptor_pvalue) ASC
        LIMIT max_results;
    END IF;
END$$

DELIMITER ;

# Example Procedure Calls
CALL GetDiseaseLigandReceptor('AD', 'Blood', 30);

# Conditions: 
# PD	Blood
# PD	BA9 prefrontal cortex
# HD	Blood
# HD	BA9 prefrontal cortex
# AD	Blood
# AD	Post mortem brain
# ALS	Cervical spinal cord
# ALS	PBMC
# ALS	Blood