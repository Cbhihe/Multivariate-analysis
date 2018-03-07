##          A Multivariate Analysis Project Proposal
####               by Cedric Bhihe, 2018.03.07


***Urban complaint calls to 311 in NYC:***

Each year, starting in 2011, between 70,000 and 120,000 calls to 311 are received by NYC's public office and logged with a slew of metadata attributes. Calls are registered with 50+ categorical, temporal and ordinal attributes, gender of caller, exact geolocation (GPS), etc. Some field attributes consist of free text pertaining to the issue resolution protocol for the calls or to the description of that issue. We do not expect to have enough time to parse natural language in order to exploit such attributes. We therefore would reduce active factors to a maximum of 25 attributes per observation. In addition, although not necessary for clustering, categorical data can be easily transformed (pre-processing) in times series over sliding windows of varying spans.

##### Year of availability:
Data is available in annual batches from 2011 to 2015 (included). 
Data sets are available at: https://nycopendata.socrata.com/Social-Services/311-Service-Requests-from-2010-to-Present/erm2-nwe9

##### Data usage and access: 
All the data is public, regulated by the terms of use of data and information available on the source web page: http://www1.nyc.gov.  The terms and conditions of use are included here and available at: http://www1.nyc.gov/home/terms-of-use.page.
The data is also available via the documented Socrata Open Data API (SODA). SODA provides programmatic access including the ability to filter, query, and aggregate data.  All communication with the API is done through HTTPS, and errors are communicated through HTTP response codes. Available response types include JSON, XML, and CSV, which are selectable by the "extension" (.json, etc.) on the API endpoint or through content-negotiation with HTTP Accepts headers.  The documentation also includes inline, runable examples.  
We do not expect to have to resort to the API to execute the MVA project.

##### Data format:
Data comes as csv, json o geojson files (one file per year) weighing about 200MB. Every dataset entry is labeled.  A rapid inspection shows that "NA" (non-assigned) values exist but in small proportion, making it possible to carry out routine imputation or suppression without incurring in issues of bias or unreliability.

###### Ideas:
In very broad terms, beyond preprocessing and validation protocols and visualization, it would be interesting to focus on:
- Clustering of registered issues and complaints
- Geo-correlative study (GPS data for each complaint is available)
- Density estimation (mechanisms which generate the data)
- Predictibility of complaints using a MVA (regression?) model, e.g.:
	- use October + November of 2013 as training set to estimate trends in december 2013, or if time permits use 2013 data to predict 2014 trends. 
	- detect micro trends that announce correlated complaints, taking seasonality as well as district location into account.

- Correlation for instance with SAT scores (extraneous to this dataset),
to explore the relationship between quantity, type and variety of complaints or disturbances and educational level of inhabitants.
