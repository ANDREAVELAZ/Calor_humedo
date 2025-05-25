pip install esgf-pyclient
pip install xarray
pip install requests

#Importing necessary libraries
from pyesgf.search import SearchConnection
import os
import pandas as pd
import requests

#Data search based on criteria

# Here we are defining criteria to be taken into account to obtain the desired dataset. 

project = ["CMIP6"] # This is to indicate that we want CMIP6 data

source_id = ["HadGEM3-GC31-MM"] # Here is to specify the models you are searching for. The list can contain more than one model.

experiment_id = ["historical", "ssp585"] # Here is to specify if you want historical data or projection (ssp245, ssp585, etc.)

variant_label = ["r1i1p1f3"] # Here is to specify the variant label

frequency = ["day"] # Here is to specify the frequency

variable = ["tas", "huss", "psl"] # Here is to specify the climate variable you are searching for

facets = ["project", "source_id", "experiment_id", "variable", "frequency", "variant_label", "latest", "replica"]

# Here we specify the node from which we want the data search to start. Setting the "distrib" parameter to "true" will expand the data search to other nodes
conn = SearchConnection('https://esgf-node.ipsl.upmc.fr/esg-search', distrib = True)

# Here we launch the search
query = conn.new_context (latest = True, # Set to "true" to only get the latest version of dataset
                          replica = False, # Set to "false" to avoid getting duplicate in the search result.
                          project = project,
                          source_id = source_id,
                       experiment_id = experiment_id,
                       variable = variable,
                       frequency = frequency,
                       variant_label = variant_label,
                       facets = facets) # This is to confirm the search criteria we have set. Must always be kept here


# Here is to obtain the total number of results the search has returned
results_count = query.hit_count 

print (f"The search has returned {results_count} results")

# Initialising an empty list that will be used to store the extracted URLs
urls = [] 

# Starting the extraction of URLs.
for i in range(results_count): # This first loop will iterate over each result

    dataset = query.search()[i] # This open a dataset 

    files_list = dataset.file_context().search() # This create a list of files contained in the opened dataset

    for file in files_list: # This loop will iterate over each file of the list to extract their URLs

        urls.append(file.download_url)

    print (f"Results {i+1} out of {results_count} processed")

# Saving the URLs in an Excel spreadsheet
df = pd.DataFrame(urls, columns = ["Links"])

df.to_excel("files_url.xlsx")

df_urls = pd.read_excel(url_file)
#df_urls_list = df_urls["Links"].list()
df_urls_list=df_urls["Links"]

# Defining a fonction for the download process

def download (url_file, output_path):
    
    # Converting the URLs table into list

    df_urls = pd.read_excel(url_file)

    df_urls_list = df_urls["Links"]

    # Starting of download

    for url in df_urls_list:

        try:
        
            file_name = str(url.rsplit("/", 1)[-1])

            file_output = output_path + file_name

            response = requests.get(url, stream = True)

            if response.status_code == 200:

                with open (file_output, 'wb') as f:

                    for chunk in response.iter_content(chunk_size=1024):

                        f.write(chunk)

                print (f"Data successfuly downloaded : {file_name}")

            else:

                print (f"Download failed for {file_name}. Status code: {response.status_code}")

        except Exception as e:

            print (f"Download failed for {file_name}. Error returned {e}")
# Data download

url_file = "/home/andreagvc/data/CMIP/files_url.xlsx" # Path to the Excel spreadsheet containing the URLs

output_path = "/home/andreagvc/data/CMIP" # Path where data will be saved

download(url_file, output_path) # Lauching the download
