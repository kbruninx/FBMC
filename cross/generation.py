# Generation levels per unit from ENTSO-e

"""
This script requires the following libraries: pandas, entsoe-py
You can install them by running the following two commands in the command line:
    pip install pandas
    pip install entsoe-py

The script can be run by typing the following command in the command line:
    python generation.py 20230329 20230331 CZ

    where the first argument is the start date, the second is the end date and the third one is the bidding zone.

For Germany, the following bidding zone codes are accepted:
    DE_50HZ, DE_AMPRION, DE_TENNET, DE_TRANSNET

    further zone codes can be found at: https://github.com/EnergieID/entsoe-py/blob/4c288c518dc38eaec18793ee4ca45321650e934c/entsoe/mappings.py#L49
"""

import sys
import pandas as pd
from entsoe import EntsoePandasClient

client = EntsoePandasClient(api_key="b18dfce9-f1e3-4d07-822f-4abd1438e602")

start_input = '20230329'
end_input = '20230331'
start = pd.Timestamp(start_input, tz='Europe/Amsterdam')
end = pd.Timestamp(end_input, tz='Europe/Amsterdam')
country_code = 'CZ'

print('Exporting generation data per unit in {0}...'.format(country_code))
print('From: {0}'.format(start))
print('To: {0}'.format(end))

df_generation = client.query_generation_per_plant(country_code, start=start, end=end)
df_generation = df_generation.tz_localize(None)

df_generation.to_excel("./generation_{}_{}_{}.xlsx".format(country_code, start_input, end_input))
print('Exported to {0}'.format("./generation_{}_{}_{}.xlsx".format(country_code, start_input, end_input)))