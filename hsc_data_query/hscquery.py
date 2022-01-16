# -*- coding: utf-8 -*-
"""
Module for accessing the Hyper Suprime-Cam Subaru Strategic Program database.

A valid account for the HSC Archive is needed to use this module.
See `HSC Online Registration
<https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new>`_.

Based on the python script developed by michitaro, NAOJ / HSC Collaboration.
[`Source <https://hsc-gitlab.mtk.nao.ac.jp/snippets/17>`_]
"""
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from builtins import input
from builtins import object
from builtins import bytes, str

import os
import json
import urllib.request, urllib.error, urllib.parse
import time
import sys
import csv
import getpass
import tempfile

from astropy import units as u
from astropy.table import Table


class QueryError(Exception):
    """
    Query error class.
    """
    pass


class HSC(object):
    """
    Main class for accessing the HSC-SSP database.

    Parameters
    ----------
    survey : str, optional
        Available surveys: 'wide', 'deep', 'udeep'.
        By default is 'wide'.
    release_version : str, optional
        For the moment, only 'pdr1' is available (Public Data Release 1)
    columns : str, optional
        List of selected columns for query results.
        See the `HSP-SSP schema <https://hsc-release.mtk.nao.ac.jp/schema/>`_
        for details. By default is 'object_id, ra, dec'.
    user : str or `None`, optional
        Account name in the HSC-SSP database. If `None`, when an ``HSC``
        object is initiated, the user can introduced the account name.
    password_env : str, optional
        The account's password can be stored in a system enviroment variable.
        By default the password is searched at ``HSCPASSW``. If this
        environment variable doesn't exist, the user is asked to introduce his
        password. Use the `password_env` option with caution, since your
        password can be easily exposed!
    """

    _version = 20181012.1
    _url = 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'
    #_url = 'https://hsc-release.mtk.nao.ac.jp/datasearch/catalog_jobs/'

    def __init__(self, survey='wide', release_version='pdr3',
                 columns='object_id, ra, dec',
                 user=None, password_env='HSCPASSW'):

        if release_version == 'pdr3':
            surveys = ['wide', 'deep', 'udeep']
        else:
            surveys = ['wide', 'dud']

        if survey not in surveys:
            error_message = 'Unknown survey: {}'
            raise ValueError(error_message.format(survey))

        user, passw = self.__login(user, password_env)

        self.credential = {'account_name': user, 'password': passw}
        self.columns = columns
        self.survey = survey
        self.release_version = release_version


    def query_region(self, coords, radius, catalog='forced'):
        """
        Returns an astropy ``Table`` object with all sources
        from catalog `catalog` within radius `radius` around
        sky position `coords`.

        Parameters
        ----------
        coords : ``SkyCoord``
            Search around this position.
        radius : ``Quantity``
            Search radius (angular units)
        catalog : str, optional
            Available options: 'forced', 'meas', 'specz', or 'random'.
            See the `HSP-SSP schema <https://hsc-release.mtk.nao.ac.jp/schema/>`_
            for details. By default is 'forced'.
        """

        catalogs = ['forced', 'meas', 'specz', 'random']

        if catalog not in catalogs:
            error_message = 'Unknown survey: {}'
            raise ValueError(error_message.format(catalog))

        table = '{}_{}.{}'.format(self.release_version, self.survey, catalog)
        data_raw = self.__cone_search(coords, radius, self.columns, table)

        data = self.__clean_fits_output(data_raw)

        return data


    def send_query(self, sql, output_format='csv',
                   output_file=None, delete_job=True):
        """
        Send an SQL query `sql`.

        If `output_file` is ``None``, a preview of the results is shown.
        Otherwise, results are saved in a file with name `output_file` and
        in the format defined by `output_format`.

        Parameters
        ----------
        sql : str
            SQL query.
        output_format : str, optional
            Available formats: 'csv', 'csv.gz', 'sqlite3', or 'fits'.
        output_file : str or ``None``
            Name of the file for storing the query results. If ``None``,
            a preview of the results is shown.
        delete_job : bool
            Delete job and results from the user space. By default is ``True``.
        """
        formats = ['csv', 'csv.gz', 'sqlite3', 'fits']

        try:
            if output_file is None:
                self.__preview(self.credential, sql, sys.stdout)

            else:
                if output_format not in formats:
                    error_message = 'Unknown output format: {}'
                    raise ValueError(error_message.format(output_format))

                job = self.__submit_job(self.credential, sql, output_format)
                self.__block_until_job_finishes(self.credential, job['id'])

                with open(output_file, 'wb') as output:
                    self.__download(self.credential, job['id'], output)

                if delete_job:
                    self.__delete_job(self.credential, job['id'])

        except urllib.error.HTTPError as error:
            if error.code == 401:
                print('invalid id or password.', file=sys.stderr)

            if error.code == 406:
                print(error.read(), file=sys.stderr)

            else:
                print(error, file=sys.stderr)

        except QueryError as error:
            print(error, file=sys.stderr)

        except KeyboardInterrupt:
            if job is not None:
                self.__job_cancel(self.credential, job['id'])

            raise


    def __login(self, user, password_env):

        if user is None:
            # user = input('HSC-SSP user: ')
            user = 'xkchen'

        password_from_envvar = os.environ.get(password_env, '')
        if password_from_envvar != '':
            passw = password_from_envvar

        else:
            # passw = getpass.getpass('password: ')
            passw = 'sJkeOASUIX3n0o5S6JIBTjpInYd0VjdIPXNENy9o'

        return user, passw


    def __cone_search(self, coords, radius,
                      columns='object_id, ra, dec',
                      table='pdr1_udeep.forced'):

        query = 'SELECT {} FROM {} WHERE coneSearch(coord, {}, {}, {})'
        query = query.format(columns, table,
                             coords.ra.deg, coords.dec.deg,
                             radius.to(u.arcsec).value)

        with tempfile.NamedTemporaryFile() as temp:
            self.send_query(query, output_format='fits', output_file=temp.name)

            temp.seek(0)
            data = Table.read(temp.name, format='fits')

        return data


    def __clean_fits_output(self, fits_table):

        # Remove isnull columns
        columns = [col for col in fits_table.colnames
                   if not col.endswith('_isnull')]

        return fits_table[columns]


    def __http_json_post(self, url, data):

        data['clientVersion'] = self._version
        post_data = json.dumps(data)
        headers = {'Content-type': 'application/json'}

        req = urllib.request.Request(url, post_data.encode(), headers)
        res = urllib.request.urlopen(req)

        return res


    def __submit_job(self, credential, sql, out_format,
                     nomail=True, skip_syntax_check=True):

        url = self._url + 'submit'
        catalog_job = {
            'sql'                     : sql,
            'out_format'              : out_format,
            'include_metainfo_to_body': True,
            'release_version'         : self.release_version,
        }

        post_data = {'credential': credential, 'catalog_job': catalog_job,
                     'nomail': nomail, 'skip_syntax_check': skip_syntax_check}

        res = self.__http_json_post(url, post_data)
        job = json.load(res)

        return job


    def __job_status(self, credential, job_id):

        url = self._url + 'status'
        post_data = {'credential': credential, 'id': job_id}

        res = self.__http_json_post(url, post_data)
        job = json.load(res)

        return job


    def __job_cancel(self, credential, job_id):

        url = self._url + 'cancel'
        post_data = {'credential': credential, 'id': job_id}

        self.__http_json_post(url, post_data)


    def __preview(self, credential, sql, out):

        url = self._url + 'preview'
        catalog_job = {
            'sql'             : sql,
            'release_version' : self.release_version,
        }

        post_data = {'credential': credential, 'catalog_job': catalog_job}
        res = self.__http_json_post(url, post_data)
        result = json.load(res)

        writer = csv.writer(out)
        for row in result['result']['rows']:
            writer.writerow(row)

        result_nrows = len(result['result']['rows'])
        if result['result']['count'] > result_nrows:
            error_message = 'only top {:d} records are displayed!'
            raise QueryError(error_message.format(result_nrows))


    def __block_until_job_finishes(self, credential, job_id):

        max_interval = 5 * 60 # sec.
        interval = 1

        while True:
            time.sleep(interval)
            job = self.__job_status(credential, job_id)

            if job['status'] == 'error':
                raise QueryError('query error: {}'.format(job['error']))

            if job['status'] == 'done':
                break

            interval *= 2
            if interval > max_interval:
                interval = max_interval


    def __download(self, credential, job_id, out):

        url = self._url + 'download'
        post_data = {'credential': credential, 'id': job_id}

        res = self.__http_json_post(url, post_data)
        buffer_size = 64 * 1<<10 # 64k

        while True:
            buf = res.read(buffer_size)
            out.write(buf)
            if len(buf) < buffer_size:
                break


    def __delete_job(self, credential, job_id):
        url = self._url + 'delete'
        post_data = {'credential': credential, 'id': job_id}

        self.__http_json_post(url, post_data)
