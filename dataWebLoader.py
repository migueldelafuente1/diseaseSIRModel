'''
Created on 31 mar. 2020

@author: Miguel
'''

import re
import os
from datetime import datetime, date
from copy import deepcopy
import json

TEMPLATE = "Highcharts.chart({}, "
GRAPHS = {'Total Cases' : "'coronavirus-cases-{}'", 
          'Daily New Cases': "'graph-cases-daily'", 
          'Active Cases': "graph-active-cases-total'",
          'Total Deaths': "'coronavirus-deaths-{}'", 
          'Daily Deaths': "'graph-deaths-daily'"}
    
class DataWebLoader(object):
    """ This class access to https://www.worldometers.info/coronavirus/ by 
    coutries and download the data for Total, Daily, Active cases and Deaths 
    from the graphs:
    
        'Total Cases'
        'Daily New Cases'
        'Active Cases'
        'Total Deaths'
        'Daily Deaths'
    
    Requires installing selenium, BeautifulSoup and Chrome webdriver
    """
    
    CHROMIUM_PATH = "C:\\Users\\Miguel\\Anaconda3\\Lib\\chromedriver.exe"
    URL = "https://www.worldometers.info/coronavirus/country/{}/"
    JSON_PATH = "data.json"
    TIMESTAMP = "DOWNLOADED TIMESTAMP"
    
    @staticmethod
    def getDataScrapper(countries='spain', save=False, download_anyway=False):
        """
        WARNING: DEPRECATED -> use getData() instead
        Args: 
        :countries <str> or list <str>, for the country names in the web
            (spain, italy, us, germany, uk ... see the web).
        :download_anyway <bool> download from the web whenever you have already 
            download the data that day or not.
        Return:
        :coutry data (dict) with the data for each table in tuples:
            (datetime.date list, int list)
        """
        # avoid unnecessary downloads
        if (not download_anyway) and DataWebLoader.__loadValuesJson(countries):
            return DataWebLoader.loadJsonData()
            
        from selenium import webdriver

        driver = webdriver.Chrome(DataWebLoader.CHROMIUM_PATH)
        
        resultDict = {}
        countries = countries if isinstance(countries, list) else list(countries)
        print(" Data Scrapper Started ",
              "--------------------------------------------------------")
        for country in countries:
            driver.get(DataWebLoader.URL.format(country))
        
            print("processing data from: ", driver.title)
            dict_values =  DataWebLoader.__loadHtmlAndProcess(driver.page_source)
            
            resultDict[country] = dict_values
        
        if save:
            DataWebLoader.save(resultDict)
            
        driver.close() 
        print(" Data Scrapper Closed ",
              "---------------------------------------------------------")
        return resultDict
    
    @staticmethod
    def getData(countries='spain', save=False, download_anyway=False):
        """
        Args: 
        :countries <str> or list <str>, for the country names in the web
            (spain, italy, us, germany, uk ... see the web).
        :download_anyway <bool> download from the web whenever you have already 
            download the data that day or not.
        Return:
        :coutry data (dict) with the data for each table in tuples:
            (datetime.date list, int list)
        """
        # avoid unnecessary downloads
        if (not download_anyway) and DataWebLoader.__loadValuesJson(countries):
            return DataWebLoader.loadJsonData()
            
        import requests
        
        resultDict = {}
        countries = countries if isinstance(countries, list) else list(countries)
        print(" Data Download Started ",
              "--------------------------------------------------------")
        for country in countries:
            response = requests.get(DataWebLoader.URL.format(country))
            if not response.ok:
                print("Error while connecting to {}: exit with status code {}/{}"
                      .format(response.url, response.status_code, response.reason))
                continue
            print("processing data from: ", response.url)
            dict_values = DataWebLoader.__loadHtmlAndProcess(response.text)
            
            resultDict[country] = dict_values
        
        if save:
            DataWebLoader.save(resultDict)
            
        print(" Data has been Downloaded ",
              "---------------------------------------------------------")
        return resultDict
    @staticmethod
    def __loadHtmlAndProcess(page_source):
        """ Load the page with selenium and get the data tags, which are in the 
        javascripts for the charts. """
        from bs4 import BeautifulSoup
        
        bs = BeautifulSoup(page_source, features='html5lib')
        
        matchs = bs.find_all('script')
        html_dict = {}
        
        for match in matchs:
            if match.attrs and match.attrs.get('type') == 'text/javascript':
                for title in GRAPHS:
                    if title in match.text:
                        html_dict[title] = str(match.text)
        
        return DataWebLoader.__processHtmlJavaScriptCharts(html_dict)
    
    @staticmethod
    def __getDateFromStr(date_str):
        date_str = date_str.replace('"', '')
        date = datetime.strptime(date_str+' 2020', '%b %d %Y')
        return date.date()
    
    @staticmethod
    def __processHtmlJavaScriptCharts(html_dict):
        """ Processing the text for each graph, transforming it in a list of 
        datetime.date and integers for the reported cases."""
        
        dict_values = {}
        for title, text in html_dict.items():
            text = text.replace('\n', '')
            text = text.replace('\\', '')
            text = text.replace(');', '')
            
            if '{' in GRAPHS[title]:
                text = text.split(TEMPLATE.format(GRAPHS[title].format('log')))[0]
                text = text.split(TEMPLATE.format(GRAPHS[title].format('linear')))[-1]
            else:
                text = text.split(TEMPLATE.format(GRAPHS[title]))[-1]
            
            x_values = re.search('categories: \[(.+?)\]', text).group(1).split(',')
            y_values = re.search('data: \[(.+?)\]', text).group(1).split(',')
            y_values = [int(y) if y!='null' else 0 for y in y_values]
            
            assert len(x_values) == len(y_values), f"length of x [{len(x_values)}] doesn't match y[{len(y_values)}]"
            x_values = list(map(DataWebLoader.__getDateFromStr, x_values))
            
            dict_values[title] = (x_values, y_values)
        
        return dict_values
    
    @staticmethod
    def save(original_resultDict):
        # convert datetime object in str date
        resultDict = deepcopy(original_resultDict)
        for country, tables in resultDict.items():
            for name, table in tables.items():
                dates, vals = table
                dates = [date.strftime("%Y %m %d") for date in dates]
                
                resultDict[country][name] = (dates, vals)
        
        resultDict[DataWebLoader.TIMESTAMP] = datetime.now().strftime("%Y %m %d : %H")
        with open(DataWebLoader.JSON_PATH, 'w') as json_file:
            json.dump(resultDict, json_file)
            
    @staticmethod
    def __loadValuesJson(countries):
        """ Return False if there is not a file with the default name or if it 
        is one, it was saved within the last 4 hours."""
        if not os.path.exists(DataWebLoader.JSON_PATH):
            return False
        
        with open(DataWebLoader.JSON_PATH) as json_file:
            resultDict = json.load(json_file)
        
        if False in [country in resultDict.keys() for country in countries]:
            return False
        if DataWebLoader.TIMESTAMP in resultDict:
            download_ts = datetime.strptime(resultDict[DataWebLoader.TIMESTAMP], 
                                            "%Y %m %d : %H")
            if (datetime.now() - download_ts).seconds < 14400:
                return True
        return False
    
    @staticmethod
    def loadJsonData(fileName=''):
        """ Load data from a previously saved json, give the fileName if it
        is not the default one. """
        
        fileName = fileName if fileName else DataWebLoader.JSON_PATH
        with open(fileName) as json_file:
            resultDict = json.load(json_file)
        
        if DataWebLoader.TIMESTAMP in resultDict:
            del resultDict[DataWebLoader.TIMESTAMP]
        # convert datetime object in string date
        for country, tables in resultDict.items():
            for name, table in tables.items():
                dates, vals = table
                dates = [datetime.strptime(d, "%Y %m %d").date() for d in dates]
                
                resultDict[country][name] = (dates, vals)
        return resultDict
        
    