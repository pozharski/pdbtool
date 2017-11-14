import os, sys, sqlite3

create_statements = [
'CREATE TABLE pdbcodes (pdbcode text PRIMARY KEY, status integer)',
 'CREATE TABLE helices (pdbcode text, helix_id integer, main_axes real, ecc1 real, ecc2 real, hlen integer, htype integer, FOREIGN KEY (pdbcode) REFERENCES pdbcodes (pdbcode) ON DELETE CASCADE ON UPDATE CASCADE, PRIMARY KEY (pdbcode, helix_id))',
 'CREATE TABLE strands (pdbcode text, sheet_id integer, strand_id integer, main_axes real, ecc1 real, ecc2 real, blen integer, btype integer, FOREIGN KEY (pdbcode) REFERENCES pdbcodes (pdbcode) ON DELETE CASCADE ON UPDATE CASCADE, PRIMARY KEY (pdbcode, sheet_id, strand_id))'
]

class sstin_dbase:
    def __init__(self, sqlite_file):
        if not os.access(sqlite_file, os.W_OK):
            conn = sqlite3.connect(sqlite_file)
            cur = conn.cursor()
            for stmt in create_statements:
                cur.execute(stmt)
            conn.commit()
            conn.close()
        self.conn = sqlite3.connect(sqlite_file)
        self.cur = self.conn.cursor()
    def fetch_processed_codes(self):
        return zip(*self.cur.execute('select pdbcode from pdbcodes where status=1').fetchall())[0]
    def filter_codes(self, codes):
        return sorted(set(codes).difference(self.fetch_processed_codes()))
    def insert_new_code(self, code):
        try:
            self.cur.execute('INSERT INTO pdbcodes (pdbcode, status) VALUES (?,0)' , tuple([code]))
            self.conn.commit()
        except sqlite3.IntegrityError:
            pass
    def insert_helix(self, code, h, t):
        try:
            self.cur.execute('INSERT INTO helices (pdbcode, helix_id, main_axes, ecc1, ecc2, hlen, htype) VALUES (?,?,?,?,?,?,?)', tuple([code,h.serNum] + t.eccent + [h.length, h.helixClass] ))
        except sqlite3.IntegrityError:
            self.cur.execute('INSERT INTO helices (pdbcode, helix_id, main_axes, ecc1, ecc2, hlen, htype) VALUES (?,?,?,?,?,?,?)', tuple([code,'a'+str(h.serNum)] + t.eccent + [h.length, h.helixClass] ))
    def insert_strand(self, code, b, t, length):
        try:
            self.cur.execute('INSERT INTO strands (pdbcode, sheet_id, strand_id, main_axes, ecc1, ecc2, blen, btype) VALUES (?,?,?,?,?,?,?,?)', tuple([code, b.sheetID, b.strand] + t.eccent + [length, b.sense]))
        except sqlite3.IntegrityError:
            self.cur.execute('INSERT INTO strands (pdbcode, sheet_id, strand_id, main_axes, ecc1, ecc2, blen, btype) VALUES (?,?,?,?,?,?,?,?)', tuple([code, b.sheetID+'1', b.strand] + t.eccent + [length, b.sense]))
    def code_lock(self, code):
        self.cur.execute('UPDATE pdbcodes SET status=1 WHERE pdbcode=?', tuple([code]))
    def commit(self):
        self.conn.commit()
    def close(self):
        self.conn.commit()
        self.conn.close()
