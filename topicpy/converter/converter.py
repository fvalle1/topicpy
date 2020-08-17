import pandas as pd
url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
df_conversion = pd.read_csv(url, sep="\t", index_col=1).dropna()

def convert_list_to_sybmols(ensgs):
	return [df_conversion.at[g,"Approved symbol"]  if g in df_conversion.index else g for g in ensgs]

def convert_list_to_ensg(symbols):
	df_inverted = df_conversion.reset_index().set_index("Approved symbol")
	return [df_inverted.at[g,"Ensembl ID(supplied by Ensembl)"]  if g in df_inverted.index else g for g in symbols]