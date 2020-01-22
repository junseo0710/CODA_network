"""
Created on 2017. 1. 18.

@author: jmjung, mjkwon, sryim
"""

import sys
sys.path.insert(0, '/data5/coda3_0/BSML')
from BSML_main import BSML_W

from datetime import datetime

import psycopg2 as pg2


def getKUTable(dictSelectCondition):
    """
    :param dictSelectCondition: a dictionary specifying selecting condition.
                                See a comment in the main function for more details.
    :return: a list of selected knowledge units
    """

    sql = "SELECT kuID, leftEntityID, leftOrganID, leftTissueID, leftCellID, association, rightEntityID, rightOrganID, rightTissueID, rightCellID, " \
          "associationContext, associationInSource, speciesID, a_referenceID, evidenceScore, level FROM KnowledgeUnit"  # attributes starting with 'A_' are of array types.

    listCondition = []  # list of sql conditions that will be used in where clause

    for key, value in dictSelectCondition.items():  # key: (tableName, columnName[, left|right])
        # If the value is an empty tuple, do not consider it
        if value == ():
            continue

        # If the value is "NOT NULL"/"NULL", then select knowledge units whose column value is NOT NULL/NULL
        elif value in ["NOT NULL", "NULL"]:
            if key[0] == "KnowledgeUnit":
                # mogrify(string, tuple) converts python variables into appropriate SQL data types.
                # mogrify() converts a string into an escape string by adding E and single quotes.
                # For example, "foo": a string. "E'foo'": an escape string.
                # To avoid this for column names, AsIs() should be used.
                if key[1] in ["associationContext", 'evidenceScore']:
                    listCondition.append(cur.mogrify("%s->'%s' IS %s", (AsIs(key[1]), AsIs(key[2]), AsIs(value))))
                else:
                    listCondition.append(cur.mogrify("%s IS %s", (AsIs(key[1]), AsIs(value))))

            elif key[0] == 'GeneProduct':
                listCondition.append(cur.mogrify("(%sEntityID IN (SELECT geneProductID FROM GeneProduct WHERE %s IS %s) OR "
                                                 "%sEntityID IN (SELECT DISTINCT geneComplexID FROM GeneComplex WHERE geneProductID IN (SELECT GeneProductID FROM GeneProduct WHERE %s IS %s)))",
                                                 (AsIs(key[2]), AsIs(key[1]), AsIs(value), AsIs(key[2]), AsIs(key[1]), AsIs(value))))

            elif key[0] == 'Reference':
                listCondition.append(cur.mogrify("A_referenceID && ARRAY(SELECT refid FROM reference WHERE %s IS %s)", (AsIs(key[1]), AsIs(value))))

        # If the value is a tuple of values that you want are searching for.
        else:
            if key[0] == "KnowledgeUnit":
                if key[1] in ["associationContext", 'evidenceScore']:
                    listCondition_sub = []
                    for entity in value:
                        if entity.isdigit():
                            sql_sub = "%s @> '{\"%s\":%s}'" % (key[1], key[2], entity)
                        else:
                            sql_sub = "%s -> '%s' ? '%s'" % (key[1], key[2], entity)
                        listCondition_sub.append(sql_sub)
                    listCondition.append(' OR '.join(listCondition_sub).encode())

                else:
                    listCondition.append(cur.mogrify("%s IN %s", (AsIs(key[1]), value)))

            elif key[0] == 'GeneProduct':
                if key[1] == 'A_modification':
                    listCondition.append(cur.mogrify("(%sEntityID IN (SELECT geneProductID FROM GeneProduct WHERE %s && %s::CHAR(1)[]) OR "
                                                     "%sEntityID IN (SELECT DISTINCT geneComplexID FROM GeneComplex WHERE geneProductID IN "
                                                     "(SELECT GeneProductID FROM GeneProduct WHERE %s && %s::CHAR(1)[])))",
                                                     (AsIs(key[2]), AsIs(key[1]), list(value), AsIs(key[2]), AsIs(key[1]), list(value))))
                else:
                    listCondition.append(cur.mogrify("(%sEntityID IN (SELECT geneProductID FROM GeneProduct WHERE %s IN %s) OR "
                                                     "%sEntityID IN (SELECT DISTINCT geneComplexID FROM GeneComplex WHERE geneProductID IN "
                                                     "(SELECT GeneProductID FROM GeneProduct WHERE %s IN %s)))",
                                                     (AsIs(key[2]), AsIs(key[1]), value, AsIs(key[2]), AsIs(key[1]), value)))

            elif key[0] == 'Reference':
                listCondition.append(cur.mogrify("A_referenceID && ARRAY(SELECT refid FROM reference WHERE %s IN %s)", (AsIs(key[1]), value)))

    if len(listCondition) > 0:
        sql += ' WHERE '+' AND '.join([cnd.decode('utf-8') for cnd in listCondition])

    # print("Query sql: ", sql)
    cur.execute(sql)

    dictKU = {}  # dictKU[kuID] = (leftEntityID, leftOrganID, leftTissueID, leftCellID, association, rightEntityID, rightOrganID, rightTissueID, rightCellID, associationContext,
                   # associationInSource, speciesID, [A_referenceID], evidenceScore, level)

    # fetchall(): list of tuples.
    # Each tuple corresponds to one row in the table and each element in a tuple corresponds to each column.
    # For example, [(col1, col2, col3, ...),  # row1
    #               (col1, col2, col3, ...)]  # row2
    listKU = cur.fetchall()
    for KU in listKU:
        dictKU[KU[0]] = list(KU[1:])

    print('The number of knowledge units: ', cur.rowcount)
    return dictKU


def getTable(tableName, PKName, setPKValue):
    """
    Return a sub-table in the form of a dictionary.
    A sub-table contains rows that appeared in the selected knowledge units at least once.
    :param tableName
    :param PKName: name of primary key of the table
    :param setPKValue: a set of primary key values that appeared in the selected knowledge units
    :return: dictTable. dictTable[primary key]=[value1, value2, ...]
    """
    # If no row appeared in the selected knowledge units, then return an empty dictionary.
    if len(setPKValue) == 0:
        return {}

    # If some row appeared in the selected knowledge units, then return the table as a dictionary.
    dictTable = {}

    if tableName == 'GeneComplex':
        cur.execute("SELECT geneComplexID, array_agg(geneProductID) FROM GeneComplex WHERE geneComplexID in %s GROUP BY geneComplexID", (tuple(setPKValue), ))
        listRow = cur.fetchall()  # List of tuples. Each tuple corresponds to a row in the table.
        for row in listRow:
            dictTable[row[0]] = row[1]

    else:
        sql = cur.mogrify("SELECT * FROM %s WHERE %s IN %s", (AsIs(tableName), AsIs(PKName), tuple(setPKValue)))
        cur.execute(sql)

        listRow = cur.fetchall()  # List of tuples. Each tuple corresponds to a row in the table.

        if tableName == 'Reference':
            for row in listRow:
                row = list(row)

                for index in [5, 6, 9]:  # version, sourceDate, exPI
                    if row[index] is None:
                        row[index] = 'NA'

                if row[10] is None:  # subEx
                    row[10] = 'NA'
                else:
                    row[10] = row[10][0]

                row[12] = row[12][0]  # BSML converter format

                dictTable[row[0]] = row[1:]

        elif tableName == 'GeneProduct':
            for row in listRow:
                dictTable[row[0]] = list(row[1:])
                if row[3] is None:
                    dictTable[row[0]][2] = 'NA'
                else:
                    dictTable[row[0]][2] = '|'.join(row[3])

                for index in [2, 4]:
                    if row[index] is None:
                        dictTable[row[0]][index - 1] = 'NA'

        else:
            for row in listRow:
                dictTable[row[0]] = list(row[1:])

    return dictTable

def completeBSMLFormat(dictKU, fileName):
    """
    Converting knowledge units into BSML format
    :param dictKU: dictKU[kuID] = (leftEntityID, leftOrganID, leftTissueID, leftCellID, association, rightEntityID, rightOrganID, rightTissueID, rightCellID, associationContext,
                                       associationInSource, speciesID, [A_referenceID], evidenceScore, level)
    :param fileName: a file to which BSML-formatted knowledge units will be written
    listKU contains information in the knowledge unit table.
    To convert knowledge units into BSML format, we also need information in other tables: GeneProduct, Reference, (GeneComplex)
    So, firstly we collect primary keys of each table that appeared in the selected knowledge units and stored them in a set.
    Secondly, for each table, we get the sub-table whose primary keys are in the set and store it in a dictionary.
    Thirdly, convert each knowledge units into BSML format by using the dictionaries.
    """
    # Set of IDs appeared in the listKU
    setRefID = set()
    setGPID = set()
    setGCID = set()

    fOut = open(fileName, 'w')

    for KU in dictKU.values():
        for entityID in [KU[0], KU[5]]:
            if entityID[:2] == 'GP':
                setGPID.add(entityID)
            elif entityID[:2] == 'GC':
                 setGCID.add(entityID)

        setRefID.update(KU[-3])

    # getTable(): get sub-table and store it in a dictionary.
    refTable = getTable('Reference', 'refID', setRefID)
    GCTable = getTable('GeneComplex', 'geneComplexID', setGCID)
    for listGP in GCTable.values():
        setGPID.update(listGP)
    GPTable = getTable('GeneProduct', 'geneProductID', setGPID)

    # Formatting each knowledge unit into a BSML format and write to a file. In BSML format, NULL is represented by 'NA'.
    for KU in dictKU.values():  # KU: tuple
        dictEntity = {}
        dictOrgan = {}
        dictTissue = {}
        dictCell = {}

        dictEntity['left'], dictOrgan['left'], dictTissue['left'], dictCell['left'], association, dictEntity['right'], dictOrgan['right'], dictTissue['right'], dictCell['right'], \
        associationContext, associationInSource, speciesID, A_referenceID, evidenceScore, level = KU

        # Formatting anatomical context
        for dictAntm in [dictOrgan, dictTissue, dictCell]:
            for key in ['left', 'right']:
                if dictAntm[key] is None:
                    dictAntm[key] = 'NA'

        # Formatting associationContext
        assoCnxt_BSML = {"Phenotype": ["NA"], "Compound": ["NA"], "Herb": ["NA"]}
        if associationContext is not None:
            assoCnxt_BSML.update(associationContext)

        # Formatting associationInSource
        if associationInSource is None:
            associationInSource = 'NA'

        # Formatting speciesID
        if speciesID is None:
            speciesID = 'NA'

        # Formatting reference
        ref_BSML = []  # list of dictionaries.
        listKey = ["1.refType", "2.name", "3.description", "4.recordID", "5.version", "6.sourceDate", "7.acqDate",
                   "8.regDate", "9.exPI", "10.subExPI","11.procPI","12.subProcPI"]

        for refID in A_referenceID:
            dictRef = {}
            for key, value in zip(listKey, refTable[refID]):
                dictRef[key] = value
            ref_BSML.append(dictRef)

        # Formatting evidenceScore
        if evidenceScore is None:
            score_BSML = 'NA'
        else:
            # score_BSML = list(evidenceScore.keys())[0] + ':' + str(list(evidenceScore.values())[0])
            score_BSML = evidenceScore

        # Formatting left event and right event
        lt_BSML = []
        rt_BSML = []
        for term_BSML, LR in zip([lt_BSML, rt_BSML], ['left', 'right']):
            entityID = dictEntity[LR]
            # entityID#geneProductType#modification#isoform,organ,tissue,cell
            if entityID[:2] == 'GP':
                entity_BSML = '%s,%s,%s,%s' % ('#'.join(GPTable[entityID]), dictOrgan[LR], dictTissue[LR], dictCell[LR])

            elif entityID[:2] == 'GC':
                listComponent = []
                for component in GCTable[entityID]:
                    component_BSML = '#'.join(GPTable[component])
                    listComponent.append(component_BSML)
                entity_BSML = '%s,%s,%s,%s' % ('&'.join(listComponent), dictOrgan[LR], dictTissue[LR], dictCell[LR])

            else:
                entity_BSML = '%s#NA#NA#NA,%s,%s,%s' % (entityID, dictOrgan[LR], dictTissue[LR], dictCell[LR])

            term_BSML.append(entity_BSML)

#        print(lt_BSML, association, rt_BSML, assoCnxt_BSML, associationInSource, speciesID, ref_BSML, score_BSML)
#         print (ref_BSML)
        if association == 'Undirected link':
            association = 'Undirected Link'
        # print(lt_BSML, association, rt_BSML, assoCnxt_BSML, associationInSource, speciesID, ref_BSML, score_BSML)
        CODAKU = BSML_W(lt_BSML, association, rt_BSML, assoCnxt_BSML, associationInSource, speciesID, ref_BSML, score_BSML)
        fOut.write(str(CODAKU) + '\t' + level + '\n')

    conn.commit()
    # fOut.write('# Process time: ' + str(datetime.now() - startTime) + '\n')
    fOut.close()


# main function
if __name__ == '__main__':
    """
    Tables used here and their attributes.
    PK: primary key.
    Attributes starting with 'A_' is of array type.

    KnowledgeUnit: kuID(PK), kugID, leftType, leftEntityID, leftOrganID, leftTissueID, leftCellID, association,
                   rightType, rightEntityID, rightOrganID, rightTissueID, rightCellID, associationContext, associationInSource, speciesID, A_referenceID, evidenceScore, level
    GeneProduct: geneProductID(PK), geneID, type, A_modification
    Reference: refID(PK), refType, name, description, recordID, version, acqDate, exPI, A_subExPI, procPI, A_subProcPI
    """

    # Connection to an existing DB. Commands using the same connection will be regarded as the same transaction.
    AsIs = pg2.extensions.AsIs
    conn = pg2.connect(database='coda3_0', user='bisler', host="heart5.kaist.ac.kr", password='bislaprom3')

    # Open a cursor to perform DB operations. The cursor is used to traverse the records.
    cur = conn.cursor()

    '''
    dictSelectCondition: a dictionary specifying conditions that will be used in the where clause of select statement.
    A Dictionary is composed of pairs of a key and a value. dictionary = {key1:value1, key2:value2:, ...}
    In dictSelectCondition, each key is a tuple composed of two or three elements: (table name, column name [, left|right])
    Each value can be in one of the four forms:
    1. An empty tuple: if the value is an empty tuple, then the column will not be considered during the knowledge unit selection.
    2. A tuple of values that you want are searching for.
       Selected knowledge units will satisfy at least one of the element in the tuple.
       In python, a tuple with only one element should include comma. For example, (1, ) is a tuple and (1) is int.
    3. A string "NOT NULL": knowledge units that have NULL for this column will not be selected.
    4. A string "NULL": knowledge units that have value for this column will not be selected.
    Each element in a value will be operated by OR operator, whereas each value will be operated by AND operator.
    For example, the following dictSelectCondition will take knowledge units that satisfy 2 conditions simultaneously:
    a. leftEntityID is 'CP00000981'
    b. rightGeneID is either 'GE16132749', or 'GE00000001'
    dictSelectCondition = {
    ("KnowledgeUnit", "leftEntityID"): ('CP00000981', ),
    ("GeneProduct", "geneID", "left"): ('GE16132749', 'GE00000001')}

    :note
    1. association = ["Positive Cause", "Positive Increase", "Negative Decrease", "Negative Cause", "Positive Decrease", "Negative Increase", "Directed Link", # directed associations
                     "Positive Correlation", "Negative Correlation", "Undirected Link"] # undirected associations
    2. reference name = ["KEGG", "EndoNet", "GO", "PhenoGO", "CTD", "BioGRID", "TRANSFAC", "ReconX", "UMLS", "RegNetwork", "PubMed", "HPA", "DiseaseConnect", # DB
                         "in vivo", "Analysis"]
    '''

    dictSelectCondition = {("KnowledgeUnit", "leftType"): (),
                           ("KnowledgeUnit", "leftEntityID"): (),
                           ("GeneProduct", "geneID", "left"): (),
                           ("GeneProduct", "type", "left"): (),
                           ("GeneProduct", "A_modification", "left"): (),
                           ("GeneProduct", "isoform", "left"): (),
                           ("KnowledgeUnit", "leftOrganID"): (),
                           ("KnowledgeUnit", "leftTissueID"): (),
                           ("KnowledgeUnit", "leftCellID"): (),
                           ("KnowledgeUnit", "association"): (),
                           ("KnowledgeUnit", "rightType"): (),
                           ("KnowledgeUnit", "rightEntityID"): ('PH00074992',),
                           ("GeneProduct", "geneID", "right"): (),
                           ("GeneProduct", "type", "right"): (),
                           ("GeneProduct", "A_modification", "right"): (),
                           ("GeneProduct", "isoform", "right"): (),
                           ("KnowledgeUnit", "rightOrganID"): (),
                           ("KnowledgeUnit", "rightTissueID"): (),
                           ("KnowledgeUnit", "rightCellID"): (),
                           ("KnowledgeUnit", "associationContext", "Herb"): (),
                           ("KnowledgeUnit", "associationContext", "Compound"): (),
                           ("KnowledgeUnit", "associationContext", "Phenotype"): (),
                           ("KnowledgeUnit", "associationInSource"): (),
                           ("KnowledgeUnit", "speciesID"): (),
                           ("Reference", "name"): (),
                           ("Reference", "refType"): (),
                           ("Reference", "recordID"): (),
                           ("Reference", "exPI"): (),
                           ("Reference", "procPI"): (),
                           ("KnowledgeUnit", 'evidenceScore', 'Manual Curation'): (),
                           ("KnowledgeUnit", 'evidenceScore', 'Literature'): (),
                           ("KnowledgeUnit", 'evidenceScore', 'P Value'): (),
                           ("KnowledgeUnit", 'evidenceScore', 'Discrete Level'): (),
                           ("KnowledgeUnit", 'evidenceScore', 'Inferred Quantity'): (),
                           ("KnowledgeUnit", "level"): ('0',)}

    startTime = datetime.now()
    print("Start: " + str(startTime))

    # select knowledge units from knowledge unit table that satisfy conditions in the dictSelectCondition.
    KUTable = getKUTable(dictSelectCondition)

    # Get more detailed information about the knowledge unit and write down the knowledge units to 'test.txt' in BSML format.
    completeBSMLFormat(KUTable, '../stomach_neoplasm_level0_right.txt')

    # Close communication with the database
    cur.close()
    conn.close()

    endTime = datetime.now()
    print("End: " + str(endTime))
    print("Process time: " + str(endTime - startTime))
