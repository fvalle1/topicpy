#  Copyright (c) 2020 fvalle
#
#  Permission is hereby granted, free of charge, to any person
#  obtaining a copy of this software and associated documentation
#  files (the "Software"), to deal in the Software without
#  restriction, including without limitation the rights to use,
#  copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following
#  conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#  OTHER DEALINGS IN THE SOFTWARE.

import pandas as pd
url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
df_conversion = pd.read_csv(url, sep="\t", index_col=1).dropna()

def convert_list_to_sybmols(ensgs: list) -> list:
	"""
	it converts a list of ENSG to gene names
    
    :param ensgs: list of ENSG

	convert_list_to_sybmols(["ENSG00000159763"])
	"""
	return [df_conversion.at[g,"Approved symbol"]  if g in df_conversion.index else g for g in ensgs]

def convert_list_to_ensg(symbols: list) -> list:
	"""
	it converts a list of gene names to ENSG identifiers

    :param symbols: list of gene symbols

	convert_list_to_sybmols(["PIP"])
	"""
	df_inverted = df_conversion.reset_index().set_index("Approved symbol")
	return [df_inverted.at[g,"Ensembl ID(supplied by Ensembl)"]  if g in df_inverted.index else g for g in symbols]
