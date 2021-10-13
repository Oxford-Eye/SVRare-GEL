import os
import yaml
import mysql.connector

# tables attributes
FILEPATH = os.path.dirname(os.path.realpath(__file__))

_default_attributes_file = os.path.join(FILEPATH, 'attributes.yml')

DEFAULT_ATTRIBUTES = {}

with open(_default_attributes_file, 'rt') as inf:
    DEFAULT_ATTRIBUTES = yaml.safe_load(inf)

# mysql connection params
_mysql_params_file = os.path.join(FILEPATH, '..', 'mysql.yml')
MYSQL_PARAMS = {}

with open(_mysql_params_file, 'rt') as inf:
    MYSQL_PARAMS = yaml.safe_load(inf)

def db_connect(params = None):
    if params is None:
        params = MYSQL_PARAMS

    con = mysql.connector.connect(**params['mysql'])
    return con
    
def create_tables(attributes = None, tables = None):
    # if tables is None, make tables for all
    result = []
    if attributes is None:
        attributes = DEFAULT_ATTRIBUTES
    if tables is None:
        tables = list(attributes.keys())

    # leave relation tables last to create, which have '_' in the name
    for table_name in sorted(tables, key=lambda x: '_' in x):
        table_values = attributes[table_name]
        # drop and create
        sql = f"DROP TABLE IF EXISTS `{table_name}`"
        result.append(sql)
        sql = f"CREATE TABLE `{table_name}` ("
        column_sql = []
        for column in table_values['columns']:
            this_sql = f"`{column['name']}` {column['type']}{column['name'] == 'id' and ' NOT NULL AUTO_INCREMENT' or ''}"
            if column['required']:
                this_sql += " NOT NULL"
            column_sql.append(this_sql)
        sql += ', '.join(column_sql)
        # add primary key
        sql += ', PRIMARY KEY (id)'
        # add foreign keys
        for foreign in table_values['foreign_keys']:
            sql += f", CONSTRAINT `fk1_{table_name}_{foreign['key']}_to_{foreign['references']['table']}_{foreign['references']['column']}` FOREIGN KEY (`{foreign['key']}`) REFERENCES `{foreign['references']['table']}` (`{foreign['references']['column']}`)"
            for cascade in foreign['cascade']:
                sql += f" ON {cascade} CASCADE"

        sql += f") ENGINE={table_values['engine']}"
        result.append(sql)
        for index in table_values['indices']:
            sql = f"CREATE INDEX index_{table_name}_{index} ON `{table_name}`({index})"
            result.append(sql)

    return result


def get_columns(conn, table):
    sql = f'SHOW COLUMNS FROM `{table}`'
    cur = conn.cursor()
    cur.execute(sql)
    columns = [i[0] for i in cur.fetchall()]
    cur.close()
    return columns

def add_record(cur, table_name, record, attributes = None):
    if attributes is None:
        attributes = DEFAULT_ATTRIBUTES
    # all required fields are present?
    required_fields = set([i['name'] for i in attributes[table_name]['columns'] if i['required']])
    if isinstance(record, dict):
        missing = required_fields - set(record.keys())
    else:
        # instance of a class
        missing = set(filter(lambda x: not hasattr(record, x), required_fields))
    if missing:
        err = f"required fields missing: {', '.join(missing)}"
        raise ValueError(err)
    # insert.
    # ignore fields not stated in attributes
    if isinstance(record, dict):
        columns = [i['name'] for i in attributes[table_name]['columns'] if i['name'] in record]
    else:
        columns = [i['name'] for i in attributes[table_name]['columns'] if hasattr(record, i['name'])]

    sql = f''' INSERT INTO `{table_name}`({','.join(columns)}) VALUES({','.join(['%s'] * len(columns))}) '''
    if isinstance(record, dict):
        cur.execute(sql, [record[h] for h in columns])
    else:
        cur.execute(sql, [getattr(record, h) for h in columns])
    return cur.lastrowid

def add_records(cur, table_name, records, attributes = None):
    # make sure all records have the same columns
    if not records:
        return
    if attributes is None:
        attributes = DEFAULT_ATTRIBUTES
    required_fields = set([i['name'] for i in attributes[table_name]['columns'] if i['required']])
    columns = None
    values = []
    for record in records:
        if isinstance(record, dict):
            this_columns = [i['name'] for i in attributes[table_name]['columns'] if i['name'] in record]
        else:
            # instance of a class
            this_columns = [i['name'] for i in attributes[table_name]['columns'] if hasattr(record, i['name'])]
        if columns is None:
            columns = this_columns
            missing = required_fields - set(columns)
            if missing:
                err = f'{table_name} required fields missing: {missing}'
                raise ValueError(err)
            if isinstance(record, dict):
                values.append([record[h] for h in columns])
            else:
                values.append([getattr(record, h) for h in columns])
            continue
        if set(this_columns) != set(columns):
            err = f'columns are not consistent across records. Please use _.add_record instead'
            raise ValueError(err)
        if isinstance(record, dict):
            values.append([record[h] for h in columns])
        else:
            values.append([getattr(record, h) for h in columns])
    sql = f''' INSERT INTO `{table_name}`({','.join(columns)}) VALUES({','.join(['%s'] * len(columns))}) '''
    cur.executemany(sql, values)

def dict_factory(cur, row):
    # https://stackoverflow.com/questions/3300464/how-can-i-get-dict-from-sqlite-query
    d = {}
    for idx, col in enumerate(cur.description):
        d[col[0]] = row[idx]
    return d


                
