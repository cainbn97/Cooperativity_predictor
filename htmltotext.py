#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Converts html file to text file
'''

from bs4 import BeautifulSoup as BS
import argparse


## parse arguments from command line
parser = argparse.ArgumentParser(description='Convert HTML to text file')
parser.add_argument('file_html', help='Give HTML file location')
parser.add_argument('file_text', help='Give text file location')
args = parser.parse_args()

file_html = args.file_html
file_text = args.file_text

with open(file_html,'r') as html:
    soup = BS(html, 'lxml')

with open(file_text,'w') as text_file:
    for line in soup.findAll('tr'):
        text = line.get_text(strip=' ', separator = '\t')
        text = text.strip(' ')
        text = text.replace('\tT','')
        text = text.replace('\tC','')
        text = text.replace('\tG','')
        text = text.replace('\tA','')
        text_file.write(text)
        text_file.write('\n')