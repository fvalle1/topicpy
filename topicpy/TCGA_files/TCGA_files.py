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

import requests as rq
import json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

cases_endpt = 'https://api.gdc.cancer.gov/cases'

fields = [
    "submitter_id",
    "case_id",
    "primary_site",
    "disease_type",
    "project.project_id",
    "diagnoses.tumor_stage",
    "diagnoses.tumor_grade",
    "diagnoses.primary_diagnosis",
    "diagnoses.classification_of_tumor",
    "annotations.classification",
    "samples.tumor_code",
    "samples.tumor_descriptor",
    "annotations.case_id",
    ]

fields = ','.join(fields)

def queryFile(idFile):
    """
    Get information for file
    
    :param str idFile: file TCGA-id
    """
    filters = {
    "op": "in",
    "content":{
        "field": "files.file_name",
        "value": [idFile]
        }
    }
    params = {
    "fields": fields,
    "filters": json.dumps(filters),
    "format": "TSV",
    "size": "1"
    }
    print("quering...%s"%idFile)
    response = rq.get(cases_endpt, params = params)
    #print(response.content.decode('utf-8'))
    r = response.content.decode("utf-8").split('\r')
    data = np.array(r[1].replace('\n','').split('\t'))
    data = data.reshape(1,len(data))
    return pd.DataFrame(data=data, columns=r[0].split('\t'), index=[0])

def makePie(df, level, c, whatToLookFor = ['disease_type']):
    fig = plt.figure(figsize=(60,15))
    ax = fig.subplots(1, len(whatToLookFor))
    for i,lookFor in enumerate(whatToLookFor):
        try:
            datatotestarr = df[lookFor].values
        except:
            datatotestarr = df[lookFor]
        utype, counts = np.unique(datatotestarr, return_counts=True)
        total = len(datatotestarr)
        try:
            labels = ['\n'.join(wrap(str(l), 20)) for l in utype]
        except:
            labels = utype
        ax[i].set_title(lookFor, fontsize=44)
        patches, texts, autotexts = ax[i].pie(counts,
                                              labels=labels,
                                              autopct=lambda p: '#:%.0f'%(p * total / 100),
                                              textprops={'fontsize':30, 'color':'white', 'wrap':True})
        for t in texts:
            t.set_fontsize(24)
            t.set_wrap(True)
            t.set_color('black')
    fig.savefig("cluster_pie_level_%d_cluster_%d.png"%(level, c))

def makeTopicPie(df, level, whatToLookFor = ['disease_type']):
    fig = plt.figure(figsize=(60,15))
    ax = fig.subplots(1, len(whatToLookFor))
    for i,lookFor in enumerate(whatToLookFor):
        datatotestarr = df[lookFor].values
        utype, counts = np.unique(datatotestarr, return_counts=True)
        total = len(datatotestarr)
        try:
            labels = ['\n'.join(wrap(str(l), 20)) for l in utype]
        except:
            labels = utype
        ax[i].set_title(lookFor, fontsize=44)
        patches, texts, autotexts = ax[i].pie(counts,
                                              labels=labels,
                                              autopct=lambda p: '#:%.0f'%(p * total / 100),
                                              textprops={'fontsize':30, 'color':'white', 'wrap':True})
        for t in texts:
            t.set_fontsize(24)
            t.set_wrap(True)
            t.set_color('black')
    fig.savefig("topic_pie_level_%d.png"%level)

def queryFiles(files):
    """
    Get infor for a list of files
    
    :param list files: list of TCGA-ids
    """
    df = pd.DataFrame(columns=fields.split(','))
    for i,f  in enumerate(files):
        df = df.append(queryFile(f), ignore_index=True, sort=True)
    #print(df.head())
    return df


def get_tcga_tissue(sample):
    """
    Get primary_site of tcga sample
    
    :param str sample: sample id
    """
    samples = pd.read_csv("/Users/filippo/Developer/tesi/results/fpkm_all/files.dat", index_col=[0], header=0)
    for fullsample in samples.index.values:
        if sample in fullsample:
            return samples.loc[fullsample,:]
