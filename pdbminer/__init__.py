import os, sys, sqlite3

create_statements = [
'CREATE TABLE pdbcodes (pdbcode text PRIMARY KEY, status integer)',
]

typarser = {
    int     : 'integer',
    float   : 'real',
    bytes   : 'text',
    str     : 'text',
}

def get_columns(item=None):
    if item:
        return [(t[0], typarser.get(type(t[1]),'text')) for t in item.__dict__.items()]
    else:
        return []

def make_create_statements(tables):
    retval = []
    for table, item, extras in tables:
        stmt = 'create table ' + table + ' ('
        stmt += ', '.join([' '.join(x) for x in [('pdbcode','text')]+get_columns(item)])
        stmt += ", FOREIGN KEY (pdbcode) REFERENCES pdbcodes (pdbcode) ON DELETE CASCADE ON UPDATE CASCADE"
        if extras:
            stmt += ", "+extras
        stmt += ')'
        retval.append(stmt)
    return retval

class general_reader(object):
    ''' Reads single value line produced and fills a class that can be then
        stored in a pdbase. '''
    def __init__(self, line=False, *args, **kwds):
        self.readline(line)
    def readline(self, line=False):
        if line:
            self.line = line
        else:
            self.line = None
        self._extra_reads()
    def _extra_reads(self):
        pass
    def report(self):
        return self.line + self._extra_report()
    def _extra_report(self):
        return ''

class pdbase(object):
    tables = []
    def __init__(self, sqlite_file):
        if not os.access(sqlite_file, os.W_OK):
            conn = sqlite3.connect(sqlite_file)
            cur = conn.cursor()
            for stmt in create_statements+make_create_statements(self.tables):
                cur.execute(stmt)
            conn.commit()
            conn.close()
        self.conn = sqlite3.connect(sqlite_file)
        self.cur = self.conn.cursor()
        self.table_item_class = dict([(x[0], x[1].__class__) for x in self.tables])
    def fetch_codes(self, status=None):
        codes = list(zip(*self.cur.execute('select pdbcode from pdbcodes'+('' if status is None else ' where status='+str(status))).fetchall()))
        if len(codes):
            return codes[0]
        else:
            return []
    def fetch_processed_codes(self):
        return self.fetch_codes(1)
    def filter_codes(self, codes):
        return sorted(set(codes).difference(self.fetch_processed_codes()))
    def insert_new_code(self, code, status=0):
        try:
            self.cur.execute('INSERT INTO pdbcodes (pdbcode, status) VALUES (?,?)' , tuple([code,status]))
            self.conn.commit()
        except sqlite3.IntegrityError:
            pass
    def code_lock(self, code):
        self.cur.execute('UPDATE pdbcodes SET status=1 WHERE pdbcode=?', tuple([code]))
    def code_unlock(self, code):
        self.cur.execute('UPDATE pdbcodes SET status=0 WHERE pdbcode=?', tuple([code]))
    def commit(self):
        self.conn.commit()
    def close(self):
        self.conn.commit()
        self.execute("VACUUM")
        self.conn.close()
    def insert_new_item(self, table, code, item):
        self.cur.execute('INSERT INTO '+table+' (pdbcode, ' + ', '.join(item.__dict__.keys())+') values ('+','.join(['?']*(1+len(item.__dict__)))+')', tuple([code]+list(item.__dict__.values())))
    def delete_items(self, table, code):
        self.cur.execute('DELETE FROM '+table+' WHERE pdbcode=?', tuple([code]))
    def delete_all(self, table):
        self.cur.execute('DELETE FROM '+table)
    def execute(self, stmt, args=None):
        if args:
            self.cur.execute(stmt, args)
        else:
            self.cur.execute(stmt)
    def executemany(self, stmt, args):
		self.cur.executemany(stmt, args)
    def get_item_number(self, table='pdbcodes'):
        return self.cur.execute('select count(*) from '+table).fetchone()[0]
    def get_items(self, table, pdbcode=None):
        if pdbcode:
            cursor = self.cur.execute('select * from '+table+' where pdbcode=?',tuple([pdbcode]))
        else:
            cursor = self.cur.execute('select * from '+table)
        keys = [t[0] for t in cursor.description][1:]
        items = []
        itemclass = self.table_item_class[table]
        for values in cursor.fetchall():
            item = itemclass()
            item.__dict__ = dict(zip(*(keys,values[1:])))
            if pdbcode:
                items.append(item)
            else:
                items.append((values[0], item))
        return items
	def get_items_codetree(self,  table):
		cursor = self.cur.execute('select * from '+table)
		keys = [t[0] for t in cursor.description][1:]
        items = {}
        itemclass = self.table_item_class[table]
        for values in cursor.fetchall():
            item = itemclass()
            item.__dict__ = dict(zip(*(keys,values[1:])))
            if values[0] in items:
				items[values[0]].append(item)
			else:
				items[values[0]] = [item]
        return items
            
    def get_present_codes(self):
        retcodes = []
        for table in self.table_item_class.keys():
            codes = self.cur.execute('select distinct(pdbcode) from '+table).fetchall()
            if len(codes):
                retcodes += list(zip(*codes))[0]
        return list(set(retcodes))

