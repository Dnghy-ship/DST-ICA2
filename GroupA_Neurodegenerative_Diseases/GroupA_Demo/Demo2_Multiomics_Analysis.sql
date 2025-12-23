/*
Overall Function: 
Get multi-omics data (RNA + Protein) for a single disease in a specific tissue
Combines RNA DEG data and protein data, filters by significance, returns top results
*/

DROP PROCEDURE IF EXISTS GetDiseaseMultiOmics;

DELIMITER $$

CREATE PROCEDURE GetDiseaseMultiOmics(
    IN disease_code VARCHAR(10),
    IN tissue_name VARCHAR(100),
    IN max_results INT
)
BEGIN
    DECLARE condition_id_val INT;
    
    
SELECT 
    condition_id
INTO condition_id_val FROM
    Condition_table
WHERE
    disease = disease_code
        AND tissue = tissue_name;
    
    IF condition_id_val IS NULL THEN
        SELECT CONCAT('No data found for ', disease_code, ' in ', tissue_name) AS message;
    ELSE
        
        SELECT 
            r.GNO,
            r.symbol AS gene_symbol,
            g.description AS gene_description,
            ROUND(r.log2FoldChange, 3) AS RNA_log2FC,
            ROUND(r.pvalue, 6) AS RNA_pvalue,
            ROUND(p.Ratio, 3) AS Protein_ratio,
            ROUND(p.p_value, 6) AS Protein_pvalue
        FROM DEG_table r
        JOIN Proteomics_table p ON r.GNO = p.GNO 
            AND r.condition_id = p.condition_id
        JOIN Gene_table g ON r.GNO = g.GNO
        WHERE r.condition_id = condition_id_val 
            AND r.pvalue < 0.1 
            AND p.p_value < 0.1
        ORDER BY (r.pvalue + p.p_value) ASC
        LIMIT max_results;
    END IF;
END$$

DELIMITER ;

# Example Procedure Calls
CALL GetDiseaseMultiOmics('PD', 'Blood', 20);

# Disease: PD
# Tissue: Blood