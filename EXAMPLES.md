# 13 Mai 2018 modified by Christopher Goyet

# Quelques exemples de commandes

Calculer le meilleur plan d'après une liste de portails:

    python2.7 makePlan.py -o -g -s --attempts 10 -d out IngressVieuxNice.txt 

Récupérer un plan déjà calculer:

    python2.7 makePlan.py -gs -d out/ -p out/vieux.pkl


# Obtenir d'autres cartes :

Avec openstreetmap:
1) Un GET pour ouvrir la session :

Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8
Accept-Encoding: gzip, deflate, br
Accept-Language: fr,fr-FR;q=0.8,en-US;q=0.5,en;q=0.3
Cache-Control: no-cache
Connection: keep-alive
Host: www.openstreetmap.org
Pragma: no-cache
Upgrade-Insecure-Requests: 1
User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:59.0) Gecko/20100101 Firefox/59.0

2)un POST pour recupérer le token:

Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8
Accept-Encoding: gzip, deflate, br
Accept-Language: fr,fr-FR;q=0.8,en-US;q=0.5,en;q=0.3
Cache-Control: no-cache
Connection: keep-alive
Content-Length: 266
Content-Type: application/x-www-form-urlencoded
Cookie: _osm_session=65e475e5beea48b346b68376faacbebf; _pk_id.1.cf09=69d965b73507dd6d.1526216147.2.1526226931.1526226339.; _pk_ref.1.cf09=%5B%22%22%2C%22%22%2C1526226339%2C%22http%3A%2F%2Fopenstreetmap.fr%2Fprojet%22%5D; _osm_welcome=hide; _pk_ses.1.cf09=*; _osm_totp_token=496279; qos_token=729826; _osm_location=7.2624%7C43.6971%7C15%7CM
Host: www.openstreetmap.org
Pragma: no-cache
Referer: https://www.openstreetmap.org/
Upgrade-Insecure-Requests: 1
User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:59.0) Gecko/20100101 Firefox/59.0

2) Un GET pour la carte : 
Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8
Accept-Encoding: gzip, deflate, br
Accept-Language: fr,fr-FR;q=0.8,en-US;q=0.5,en;q=0.3
Cache-Control: no-cache
Connection: keep-alive
Cookie: _osm_totp_token=496279; qos_token=729826
Host: render.openstreetmap.org
Pragma: no-cache
Referer: https://www.openstreetmap.org/
Upgrade-Insecure-Requests: 1
User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:59.0) Gecko/20100101 Firefox/59.0

https://render.openstreetmap.org/cgi-bin/export?bbox=7.2686147689819345,43.69412811532704,7.281103134155274,43.69955814990594&scale=10000&format=png