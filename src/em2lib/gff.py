"""
Some utilities for GFF data manipulation in gffpandas DataFrames.
"""
#  CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21
#  Frédéric PLEWNIAK, CNRS/Université de Strasbourg UMR7156 - GMGM
#
import gffpandas.gffpandas as gffpd


def select_features(gff, references=None, ftype=None):
    """
    Feature selection from a GFF file according to reference sequence and/or type

    :param gff: the GFF file name
    :param references: a reference name or a list thereof
    :param ftype: the type of features to select
    :return: a DataFrame containing the selected features
    """
    if ftype is not None and not isinstance(ftype, list):
        ftype = [ftype]
    # ensure references is None or a list
    if references is not None and not isinstance(references, list):
        references = [references]
    gff = gffpd.read_gff3(gff)
    # remove lines which are not features (Fasta sequence, etc.)
    gff.df.dropna(axis=0, how='any', subset=['start', 'end'], inplace=True)
    # then keep only features of type specified by ftype if ftype is not None
    gff_df = gff.attributes_to_columns() if ftype is None else gff.filter_feature_of_type(ftype).attributes_to_columns()
    if references is not None:
        # keep only features that are in references
        gff_df = gff_df.loc[gff_df['seq_id'].isin(references)]
    # make sure start and end are integers
    gff_df = gff_df.astype({'start': int, 'end': int})
    return gff_df
