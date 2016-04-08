def setup(app):
    from datetime import datetime
    from random import randint
    import json
    import urllib2
    u = ('http://gdata.youtube.com/feeds/api/playlists/' +
         'UUJR-f8bqXzPv5aUlu1-Ba0g?alt=json')
    d = json.loads(urllib2.urlopen(u).read())
    n = len(d['feed']['entry'])
    i = randint(0, n - 1)
    # avoid non-ascii names of people ...
    title = d['feed']['entry'][i]['title']['$t'].split(':')[-1]
    event = d['feed']['entry'][i]['content']['$t'].split('\n')[-1]
    recorded = datetime.strptime(d['feed']['entry'][i]['yt$recorded']['$t'],
                                 '%Y-%m-%d').strftime('%d %B %Y')
    href = d['feed']['entry'][i]['link'][0]['href']

    f = open('ytp.txt', 'w')
    f.write('.. raw:: html\n')
    f.write('\n')
    f.write('       <p></p>\n')
    f.write('        <p>' + title.encode('ascii', 'ignore') + '</p>\n')
    f.write('        <p>Presented ' + recorded + ' at <a href="' +
            event.encode('ascii', 'ignore') + '">' +
            event.encode('ascii', 'ignore') + '</a>' + '</p>\n')
    f.write('       <iframe width="426" height="240" src="' +
            href.replace('watch?v=', 'embed/').replace(
                '&feature=youtube_gdata', '').replace('http:', 'https:') +
            '" frameborder="0"></iframe>\n')
    f.write('       <p></p>')
    f.close()
