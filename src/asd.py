
import pandas as pd
import bs4
import requests

heads = {'user-agent': 'Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36'}
url = "https://www.google.com/search?q=%22Barn+Owl%22&tbm=nws&tbs=qdr:y"

Links = []
headlines = []
headers = []
req = requests.get(url, headers=heads)
soup = bs4.BeautifulSoup(req.text, "html.parser")
Link  = [a["href"] for a in soup.find_all("a", class_="_HId")]
headers += len(Link) * ["Header"]
len_before=len(Link)
Link += [a["href"] for a in soup.find_all("a", class_="_sQb")]
Links.append(Link)
headers += (len(Link)-len_before) * ["Empty"]
headline =  [a.get_text() for a in soup.find_all("a", class_="_HId")]
headline += [a.get_text() for a in soup.find_all("a", class_="_sQb")]
headlines.append(headline)
df = pd.DataFrame({'headlines' : headlines[0],'Header or Not' : headers,'URLs' : Links[0],})
print (df)