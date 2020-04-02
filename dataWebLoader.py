'''
Created on 31 mar. 2020

@author: Miguel
'''

import re
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
    
    @staticmethod
    def getData(countries='spain', save=False):
        """
        Args: 
        :countries <str> or list <str>, for the country names in the web
            (spain, italy, us, germany, uk ... see the web).
        Return:
        :coutry data (dict) with the data for each table in tuples:
            (datetime.date list, int list)
        """
        from selenium import webdriver

        driver = webdriver.Chrome(DataWebLoader.CHROMIUM_PATH)
        
        resultDict = {}
        countries = countries if isinstance(countries, list) else list(countries)
        print(" Data Scrapper Started --------------------------------------------------------")
        for country in countries:
            driver.get(DataWebLoader.URL.format(country))
        
            print("processing data from: ", driver.title)
            dict_values =  DataWebLoader.__loadHtmlAndProcess(driver.page_source)
            
            resultDict[country] = dict_values
        
        if save:
            DataWebLoader.save(resultDict)
            
        driver.close() 
        print(" Data Scrapper Closed ---------------------------------------------------------")
        return resultDict
    
    @staticmethod
    def __loadHtmlAndProcess(page_source):
        """ Load the page with selenium and get the data tags, which are in the 
        javascripts for the charts. """
        from bs4 import BeautifulSoup
        
        html_source = page_source
        bs = BeautifulSoup(html_source, features='html5lib')
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
                
        with open(DataWebLoader.JSON_PATH, 'w') as json_file:
            json.dump(resultDict, json_file)
            
    @staticmethod
    def loadJsonData(fileName=''):
        """ Load data from a previously saved json, give the fileName if it
        is not the default one. """
        
        fileName = fileName if fileName else DataWebLoader.JSON_PATH
        with open(fileName) as json_file:
            resultDict = json.load(json_file)
            
        # convert datetime object in str date
        for country, tables in resultDict.items():
            for name, table in tables.items():
                dates, vals = table
                dates = [datetime.strptime(d, "%Y %m %d").date() for d in dates]
                
                resultDict[country][name] = (dates, vals)
        return resultDict
        
    