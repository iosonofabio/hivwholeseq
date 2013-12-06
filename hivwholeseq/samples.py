# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV samples from patients.

            Note: it is possible to read the Excel XLSX table with xlrd. It also
            works ok, but it's painful and error-prone code, and it assumes that
            table is up-to-date. Better write down the machine-relevant info
            manually down here (until otherwise specified).
'''
# Globals
# Date is YYYY-M-D
samples = {'VK04-3106': {'run': 37, 'date': '2004-7-9'},
           '08HR-0235': {'run': 37, 'date': '2008-2-7'},
           'VK07-4778': {'run': 37, 'date': '2007-7-12'},
           'VK03-4298': {'run': 37, 'date': '2003-9-30'},
           'VK09-7738': {'run': 37, 'date': '2009-9-24'},
           'VK08-8014': {'run': 37, 'date': '2008-10-28'},
           'VK08-6634': {'run': 37, 'date': '2008-9-11'},
           'VK07-8262': {'run': 37, 'date': '2007-11-27'},
           'VK08-2987': {'run': 37, 'date': '2008-4-21'},
}



# Functions
def date_to_integer(date):
    '''Convert a date in the format YYYY-M-D into an integer'''
    import datetime as dt
    return dt.date.toordinal(dt.datetime(*map(int, date.split('-'))))
