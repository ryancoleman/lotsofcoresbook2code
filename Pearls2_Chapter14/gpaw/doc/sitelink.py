from docutils import nodes

# Note: set_classes might have something to do with class=".." and css
# http://docutils.sourceforge.net/docs/howto/rst-roles.html

def make_site_linker(baseurl):
    def role(role, rawtext, text, lineno, inliner, options={}, content=[]):
        if text.endswith('>'):
            startindex = text.find('<')
            ref = text[startindex + 1:-1]
            label = text[0:startindex].strip()
        else:
            label = baseurl + text
            ref = text
        node = nodes.reference(rawtext, label, refuri=baseurl + ref,
                               **options)
        return [node], []
    return role

def setup(app):
    app.add_role('ase', make_site_linker('https://wiki.fysik.dtu.dk/ase/'))
    app.add_role('wiki', make_site_linker('https://wiki.fysik.dtu.dk/'))
