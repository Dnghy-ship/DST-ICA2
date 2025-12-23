/*
Overall Function:
This script has 4 stored procedures for bioinformatics analysis of common features 
across diseases in the same tissue:
1. FindCommonDEGsSameTissueTwoDiseases: Common Differentially Expressed Genes (DEGs) for 2 diseases
2. FindCommonDEGsSameTissueThreeDiseases: Common DEGs for 3 diseases
3. FindCommonKEGGPathwaysSameTissueTwoDiseases: Common KEGG pathways for 2 diseases
4. FindCommonKEGGPathwaysSameTissueThreeDiseases: Common KEGG pathways for 3 diseases
Last part: Example calls to test the procedures

Notice: The Error Code could be ignored while the demo has output in the last result part.
*/

# Procedure 1: Common DEGs (2 Diseases, Same Tissue)
DROP PROCEDURE IF EXISTS FindCommonDEGsSameTissueTwoDiseases;

DELIMITER $$

CREATE PROCEDURE FindCommonDEGsSameTissueTwoDiseases(
    IN disease1 VARCHAR(10),
    IN disease2 VARCHAR(10),
    IN tissue_name VARCHAR(100),  
    IN max_results INT
)
BEGIN
    DECLARE condition_id1 INT;
    DECLARE condition_id2 INT;
    

SELECT 
    condition_id
INTO condition_id1 FROM
    Condition_table
WHERE
    disease = disease1
        AND tissue = tissue_name;
    
SELECT 
    condition_id
INTO condition_id2 FROM
    Condition_table
WHERE
    disease = disease2
        AND tissue = tissue_name;
    

    IF condition_id1 IS NULL THEN
        SELECT CONCAT('ERROR: No ', tissue_name, ' tissue data for disease "', disease1, '"') AS error_message;
    ELSEIF condition_id2 IS NULL THEN
        SELECT CONCAT('ERROR: No ', tissue_name, ' tissue data for disease "', disease2, '"') AS error_message;
    ELSE
        
        SELECT 
            CONCAT('Comparison in ', tissue_name, ' tissue: ', disease1, ' vs ', disease2) AS analysis_title,
            CONCAT('Condition IDs: ', disease1, '=', condition_id1, ', ', disease2, '=', condition_id2) AS condition_info,
            (SELECT COUNT(*) FROM DEG_table WHERE condition_id = condition_id1 AND pvalue < 0.05) AS disease1_sig_genes,
            (SELECT COUNT(*) FROM DEG_table WHERE condition_id = condition_id2 AND pvalue < 0.05) AS disease2_sig_genes;
        
        
        SET @query_sql = CONCAT(
            'SELECT ',
            '    d1.GNO, ',
            '    d1.symbol AS gene_symbol, ',
            '    g.description AS gene_description, ',
            '    CASE ',
            '        WHEN g.is_secreted = 1 AND g.is_receptor = 1 THEN ''Secreted/Receptor'' ',
            '        WHEN g.is_secreted = 1 THEN ''Secreted Protein'' ',
            '        WHEN g.is_receptor = 1 THEN ''Receptor'' ',
            '        ELSE ''Other'' ',
            '    END AS gene_type, ',
            '    ROUND(d1.log2FoldChange, 3) AS `', disease1, '_log2FC`, ',
            '    ROUND(d1.pvalue, 6) AS `', disease1, '_pvalue`, ',
            '    ROUND(d2.log2FoldChange, 3) AS `', disease2, '_log2FC`, ',
            '    ROUND(d2.pvalue, 6) AS `', disease2, '_pvalue`, ',
            '    g.chromosome, ',
            '    g.Entry AS uniprot_entry, ',
            '    g.Protein_names AS protein_names ',
            'FROM DEG_table d1 ',
            'JOIN DEG_table d2 ON d1.GNO = d2.GNO ',
            'JOIN Gene_table g ON d1.GNO = g.GNO ',
            'WHERE d1.condition_id = ', condition_id1, ' ',
            '  AND d2.condition_id = ', condition_id2, ' ',
            '  AND d1.pvalue < 0.1 ',
            '  AND d2.pvalue < 0.1 ',
            'ORDER BY (d1.pvalue + d2.pvalue) ASC '
        );
        
        IF max_results > 0 THEN
            SET @query_sql = CONCAT(@query_sql, 'LIMIT ', max_results);
        END IF;
        
        PREPARE stmt FROM @query_sql;
        EXECUTE stmt;
        DEALLOCATE PREPARE stmt;
    END IF;
END$$

DELIMITER ;

# Procedure 2: Common DEGs (3 Diseases, Same Tissue)
DROP PROCEDURE IF EXISTS FindCommonDEGsSameTissueThreeDiseases;

DELIMITER $$

CREATE PROCEDURE FindCommonDEGsSameTissueThreeDiseases(
    IN disease1 VARCHAR(10),
    IN disease2 VARCHAR(10),
    IN disease3 VARCHAR(10),
    IN tissue_name VARCHAR(100),
    IN max_results INT
)
BEGIN
    DECLARE condition_id1 INT;
    DECLARE condition_id2 INT;
    DECLARE condition_id3 INT;
    DECLARE error_msg VARCHAR(500);
    
    
    SELECT condition_id INTO condition_id1 FROM Condition_table WHERE disease = disease1 AND tissue = tissue_name;
    SELECT condition_id INTO condition_id2 FROM Condition_table WHERE disease = disease2 AND tissue = tissue_name;
    SELECT condition_id INTO condition_id3 FROM Condition_table WHERE disease = disease3 AND tissue = tissue_name;
    
    
    SET error_msg = '';
    IF condition_id1 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' data for ', disease1, '. '); END IF;
    IF condition_id2 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' data for ', disease2, '. '); END IF;
    IF condition_id3 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' data for ', disease3, '. '); END IF;
    
    IF error_msg != '' THEN
        SELECT CONCAT('ERROR: ', error_msg) AS error_message;
    ELSE
        
        SELECT 
            CONCAT('Three-disease comparison in ', tissue_name, ' tissue') AS analysis_title,
            CONCAT(disease1, ' (ID:', condition_id1, ')') AS disease1_info,
            CONCAT(disease2, ' (ID:', condition_id2, ')') AS disease2_info,
            CONCAT(disease3, ' (ID:', condition_id3, ')') AS disease3_info;
        
        
        SET @query_sql = CONCAT(
            'SELECT ',
            '    d1.GNO, ',
            '    d1.symbol AS gene_symbol, ',
            '    g.description AS gene_description, ',
            '    CASE ',
            '        WHEN g.is_secreted = 1 AND g.is_receptor = 1 THEN ''Secreted/Receptor'' ',
            '        WHEN g.is_secreted = 1 THEN ''Secreted Protein'' ',
            '        WHEN g.is_receptor = 1 THEN ''Receptor'' ',
            '        ELSE ''Other'' ',
            '    END AS gene_type, ',
            '    ROUND(d1.log2FoldChange, 3) AS `', disease1, '_log2FC`, ',
            '    ROUND(d1.pvalue, 6) AS `', disease1, '_pvalue`, ',
            '    ROUND(d2.log2FoldChange, 3) AS `', disease2, '_log2FC`, ',
            '    ROUND(d2.pvalue, 6) AS `', disease2, '_pvalue`, ',
            '    ROUND(d3.log2FoldChange, 3) AS `', disease3, '_log2FC`, ',
            '    ROUND(d3.pvalue, 6) AS `', disease3, '_pvalue`, ',
            '    g.chromosome, ',
            '    g.Entry AS uniprot_entry ',
            'FROM DEG_table d1 ',
            'JOIN DEG_table d2 ON d1.GNO = d2.GNO ',
            'JOIN DEG_table d3 ON d1.GNO = d3.GNO ',
            'JOIN Gene_table g ON d1.GNO = g.GNO ',
            'WHERE d1.condition_id = ', condition_id1, ' ',
            '  AND d2.condition_id = ', condition_id2, ' ',
            '  AND d3.condition_id = ', condition_id3, ' ',
            '  AND d1.pvalue < 0.1 ',
            '  AND d2.pvalue < 0.1 ',
            '  AND d3.pvalue < 0.1 ',
            'ORDER BY (d1.pvalue + d2.pvalue + d3.pvalue) ASC '
        );
        
        IF max_results > 0 THEN
            SET @query_sql = CONCAT(@query_sql, 'LIMIT ', max_results);
        END IF;
        
        PREPARE stmt FROM @query_sql;
        EXECUTE stmt;
        DEALLOCATE PREPARE stmt;
    END IF;
END$$

DELIMITER ;

# Procedure 3: Common KEGG Pathways (2 Diseases, Same Tissue)
DROP PROCEDURE IF EXISTS FindCommonKEGGPathwaysSameTissueTwoDiseases;

DELIMITER $$

CREATE PROCEDURE FindCommonKEGGPathwaysSameTissueTwoDiseases(
    IN disease1 VARCHAR(10),
    IN disease2 VARCHAR(10),
    IN tissue_name VARCHAR(100),
    IN max_results INT
)
BEGIN
    DECLARE condition_id1 INT;
    DECLARE condition_id2 INT;
    
    SELECT condition_id INTO condition_id1 FROM Condition_table WHERE disease = disease1 AND tissue = tissue_name;
    SELECT condition_id INTO condition_id2 FROM Condition_table WHERE disease = disease2 AND tissue = tissue_name;
    
    IF condition_id1 IS NULL THEN
        SELECT CONCAT('ERROR: No ', tissue_name, ' tissue data for ', disease1) AS error_message;
    ELSEIF condition_id2 IS NULL THEN
        SELECT CONCAT('ERROR: No ', tissue_name, ' tissue data for ', disease2) AS error_message;
    ELSE
        
        SELECT 
            CONCAT('KEGG Pathway comparison in ', tissue_name, ' tissue') AS analysis_title,
            CONCAT(disease1, ' vs ', disease2) AS comparison,
            CONCAT('Condition IDs: ', condition_id1, ' and ', condition_id2) AS condition_info;
        
        SET @query_sql = CONCAT(
            'SELECT ',
            '    k1.PathwayID, ',
            '    k1.Description AS pathway_description, ',
            '    k1.category AS pathway_category, ',
            '    k1.subcategory AS pathway_subcategory, ',
            '    k1.GeneRatio AS `', disease1, '_gene_ratio`, ',
            '    k2.GeneRatio AS `', disease2, '_gene_ratio`, ',
            '    ROUND(k1.RichFactor, 3) AS `', disease1, '_rich_factor`, ',
            '    ROUND(k2.RichFactor, 3) AS `', disease2, '_rich_factor`, ',
            '    ROUND(k1.pvalue, 6) AS `', disease1, '_pvalue`, ',
            '    ROUND(k2.pvalue, 6) AS `', disease2, '_pvalue`, ',
            '    k1.Count AS `', disease1, '_gene_count`, ',
            '    k2.Count AS `', disease2, '_gene_count`, ',
            '    CASE ',
            '        WHEN k1.pvalue < 0.001 AND k2.pvalue < 0.001 THEN ''Both Highly Significant'' ',
            '        WHEN k1.pvalue < 0.01 AND k2.pvalue < 0.01 THEN ''Both Very Significant'' ',
            '        WHEN k1.pvalue < 0.05 AND k2.pvalue < 0.05 THEN ''Both Significant'' ',
            '        ELSE ''Partially Significant'' ',
            '    END AS significance_status ',
            'FROM KEGG_analysis k1 ',
            'JOIN KEGG_analysis k2 ON k1.PathwayID = k2.PathwayID ',
            'WHERE k1.condition_id = ', condition_id1, ' ',
            '  AND k2.condition_id = ', condition_id2, ' ',
            '  AND k1.pvalue < 0.05 ',
            '  AND k2.pvalue < 0.05 ',
            'ORDER BY (k1.pvalue + k2.pvalue) ASC '
        );
        
        IF max_results > 0 THEN
            SET @query_sql = CONCAT(@query_sql, 'LIMIT ', max_results);
        END IF;
        
        PREPARE stmt FROM @query_sql;
        EXECUTE stmt;
        DEALLOCATE PREPARE stmt;
    END IF;
END$$

DELIMITER ;

# Procedure 4: Common KEGG Pathways (3 Diseases, Same Tissue)
DROP PROCEDURE IF EXISTS FindCommonKEGGPathwaysSameTissueThreeDiseases;

DELIMITER $$

CREATE PROCEDURE FindCommonKEGGPathwaysSameTissueThreeDiseases(
    IN disease1 VARCHAR(10),
    IN disease2 VARCHAR(10),
    IN disease3 VARCHAR(10),
    IN tissue_name VARCHAR(100),
    IN max_results INT
)
BEGIN
    DECLARE condition_id1 INT;
    DECLARE condition_id2 INT;
    DECLARE condition_id3 INT;
    DECLARE error_msg VARCHAR(500);
    
    
    SELECT condition_id INTO condition_id1 FROM Condition_table WHERE disease = disease1 AND tissue = tissue_name;
    SELECT condition_id INTO condition_id2 FROM Condition_table WHERE disease = disease2 AND tissue = tissue_name;
    SELECT condition_id INTO condition_id3 FROM Condition_table WHERE disease = disease3 AND tissue = tissue_name;
    
    
    SET error_msg = '';
    IF condition_id1 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' tissue data for ', disease1, '. '); END IF;
    IF condition_id2 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' tissue data for ', disease2, '. '); END IF;
    IF condition_id3 IS NULL THEN SET error_msg = CONCAT(error_msg, 'No ', tissue_name, ' tissue data for ', disease3, '. '); END IF;
    
    IF error_msg != '' THEN
        SELECT CONCAT('ERROR: ', error_msg) AS error_message;
    ELSE
        
        SELECT 
            CONCAT('KEGG Pathway comparison in ', tissue_name, ' tissue') AS analysis_title,
            CONCAT('Comparing: ', disease1, ', ', disease2, ', ', disease3) AS diseases,
            CONCAT('Condition IDs: ', condition_id1, ', ', condition_id2, ', ', condition_id3) AS condition_info;
        
        
        SET @query_sql = CONCAT(
            'SELECT ',
            '    k1.PathwayID, ',
            '    k1.Description AS pathway_description, ',
            '    k1.category AS pathway_category, ',
            '    k1.subcategory AS pathway_subcategory, ',
            '    k1.GeneRatio AS `', disease1, '_gene_ratio`, ',
            '    k2.GeneRatio AS `', disease2, '_gene_ratio`, ',
            '    k3.GeneRatio AS `', disease3, '_gene_ratio`, ',
            '    ROUND(k1.RichFactor, 3) AS `', disease1, '_rich_factor`, ',
            '    ROUND(k2.RichFactor, 3) AS `', disease2, '_rich_factor`, ',
            '    ROUND(k3.RichFactor, 3) AS `', disease3, '_rich_factor`, ',
            '    ROUND(k1.pvalue, 6) AS `', disease1, '_pvalue`, ',
            '    ROUND(k2.pvalue, 6) AS `', disease2, '_pvalue`, ',
            '    ROUND(k3.pvalue, 6) AS `', disease3, '_pvalue`, ',
            '    k1.Count AS `', disease1, '_gene_count`, ',
            '    k2.Count AS `', disease2, '_gene_count`, ',
            '    k3.Count AS `', disease3, '_gene_count`, ',
            '    CASE ',
            '        WHEN k1.pvalue < 0.001 AND k2.pvalue < 0.001 AND k3.pvalue < 0.001 THEN ''All Highly Significant'' ',
            '        WHEN k1.pvalue < 0.01 AND k2.pvalue < 0.01 AND k3.pvalue < 0.01 THEN ''All Very Significant'' ',
            '        WHEN k1.pvalue < 0.05 AND k2.pvalue < 0.05 AND k3.pvalue < 0.05 THEN ''All Significant'' ',
            '        ELSE ''Partially Significant'' ',
            '    END AS significance_status ',
            'FROM KEGG_analysis k1 ',
            'JOIN KEGG_analysis k2 ON k1.PathwayID = k2.PathwayID ',
            'JOIN KEGG_analysis k3 ON k1.PathwayID = k3.PathwayID ',
            'WHERE k1.condition_id = ', condition_id1, ' ',
            '  AND k2.condition_id = ', condition_id2, ' ',
            '  AND k3.condition_id = ', condition_id3, ' ',
            '  AND k1.pvalue < 0.05 ',
            '  AND k2.pvalue < 0.05 ',
            '  AND k3.pvalue < 0.05 ',
            'ORDER BY (k1.pvalue + k2.pvalue + k3.pvalue) ASC '
        );
        
        IF max_results > 0 THEN
            SET @query_sql = CONCAT(@query_sql, 'LIMIT ', max_results);
        END IF;
        
        PREPARE stmt FROM @query_sql;
        EXECUTE stmt;
        DEALLOCATE PREPARE stmt;
    END IF;
END$$

DELIMITER ;

# Example Procedure Calls
CALL FindCommonDEGsSameTissueTwoDiseases('AD', 'PD', 'Blood', 30);

CALL FindCommonDEGsSameTissueThreeDiseases('AD', 'PD', 'ALS', 'Blood', 20);

CALL FindCommonKEGGPathwaysSameTissueTwoDiseases('AD', 'PD', 'Blood', 15);

CALL FindCommonKEGGPathwaysSameTissueThreeDiseases('AD', 'PD', 'HD','Blood', 15);

# Diseases: AD, ALS, PD, HD
# Tissue: Blood