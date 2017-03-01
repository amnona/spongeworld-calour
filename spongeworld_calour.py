import requests
import webbrowser
import os
import operator
from logging import getLogger

import scipy.stats

from calour.util import get_config_value
from calour.database import Database

logger = getLogger(__name__)


class SpongeWorld(Database):
    def __init__(self):
        super().__init__(database_name='SpongeWorld', methods=['get'])

        # Web address of the bact server
        self.dburl = self._get_db_address()
        self.username = get_config_value('username', section='dbBact')
        self.password = get_config_value('password', section='dbBact')
        self.web_interface = self.dburl

    def _get_db_address(self):
        '''
        Get the database address based on the environment variable SCDB_WEBSITE_TYPE
        (use export SCDB_WEBSITE_TYPE="local" / "main"(default) / "develop")

        Returns
        -------
        server_address : str
            the supercooldb server web address based on the env. variable
        '''
        if 'SPONGEWORLD_SERVER_TYPE' in os.environ:
            servertype = os.environ['SPONGEWORLD_SERVER_TYPE'].lower()
            if servertype == 'local':
                logger.debug('sponge servertype is local')
                server_address = 'http://127.0.0.1:5000'
            elif servertype == 'main':
                logger.debug('servertype is main')
                server_address = 'http://amnonim.webfactional.com/spongeworld/'
            elif servertype == 'develop':
                logger.debug('servertype is develop')
                server_address = 'http://amnonim.webfactional.com/spongeworld_develop'
            else:
                raise ValueError('unknown server type %s in SPONGEWORLD_SERVER_TYPE' % servertype)
        else:
            server_address = 'http://amnonim.webfactional.com/spongeworld'
            logger.debug('using default server main (use env. variable SPONGEWORLD_SERVER_TYPE to set)')
        return server_address

    def _post(self, api, rdata):
        '''POST a request to spongeworld and show error

        Parameters
        ----------
        api : str
            the REST API address to post the request to
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        res = requests.post(self.dburl + '/' + api, json=rdata)
        if res.status_code != 200:
            logger.warn('REST error %s enountered when accessing SpongeWorld %s: %s' % (res.reason, api, res.content))
        return res

    def _get(self, api, rdata):
        '''GET a request to SpongeWorld

        Parameters
        ----------
        api : str
            the REST API address to post the request to
        rdata : dict
            parameters to pass to the dbBact REST API

        Returns
        -------
        res : request
            the result of the request
        '''
        res = requests.get(self.dburl + '/' + api, json=rdata)
        if res.status_code != 200:
            logger.warn('REST error %s enountered when accessing SpongeWorld %s: %s' % (res.reason, api, res.content))
        return res

    def get_annotation_string(self, info, pval=0.1):
        '''Get nice string summaries of annotations

        Parameters
        ----------
        info : dict (see get_sequence_annotations)
            'total_samples' : int
                the total amount of samples in the database
            'total_observed' : int
                the total number of samples where the sequence is present
            'info' : dict of {field(str): information(dict)}
                the frequency of the sequence in each field.
                information is a dict of {value(str): distribution(dict)}
                distribution contains the following key/values:
                    'total_samples': int
                        the total number of samples having this value
                    'observed_samples': int
                        the number of samples with this value which have the sequence present in them

        Returns
        -------
        desc : list of str
            a short summary of each annotation, sorted by importance
        '''
        keep = []
        total_observed = info['total_observed']
        if total_observed == 0:
            return []
        if total_observed is None:
            logger.debug('sequence %s not found in database')
            return []
        total_samples = info['total_samples']
        null_pv = 1 - (total_observed / total_samples)
        for cfield in info['info'].keys():
            for cval, cdist in info['info'][cfield].items():
                observed_val_samples = cdist['observed_samples']
                total_val_samples = cdist['total_samples']
                cfrac = observed_val_samples / total_val_samples
                cpval = scipy.stats.binom.cdf(total_val_samples - observed_val_samples, total_val_samples, null_pv)
                if cpval <= pval:
                    cdesc = '%s:%s (%d/%d)' % (cfield, cval, observed_val_samples, total_val_samples)
                    keep.append([cdesc, cfrac, cpval])
        logger.debug('found %d significant annotations' % len(keep))

        # sort first by p-value and then by fraction (so fraction is more important)
        keep = sorted(keep, key=operator.itemgetter(2), reverse=False)
        keep = sorted(keep, key=operator.itemgetter(1), reverse=True)
        desc = [ckeep[0] for ckeep in keep]
        desc = ['Found in %f samples (%d / %d)' % (total_observed / total_samples, total_observed, total_samples)] + desc
        return desc

    def get_seq_annotations(self, sequence):
        '''Get the annotations for a sequence

        Parameters
        ----------
        sequence : str
            The DNA sequence to get the annotations for

        Returns
        -------
        annotations : list of list of (annotation dict,list of [Type,Value] of annotation details)
            See dbBact sequences/get_annotations REST API documentation
        term_info : dict of {term: info}
            where key (str) is the ontology term, info is a dict of details containing:
                'total_annotations' : total annotations having this term in the database
                'total_sequences' : number of annotations with this term for the sequence
        '''
        sequence = sequence.upper()
        rdata = {}
        rdata['sequence'] = sequence
        res = self._get('sequence/info', rdata=rdata)
        if res.status_code != 200:
            return []
        info = res.json()
        desc = self.get_annotation_string(info)
        return desc

    def get_seq_annotation_strings(self, sequence):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        sequence : str
            the DNA sequence to query the annotation strings about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails,annotationsummary)
            a list of:
                annotationdetails : dict
                    'seqid' : str, the sequence annotated
                    'annotationtype : str
                    ...
                annotationsummary : str
                    a short summary of the annotation
        '''
        shortdesc = []
        annotations = self.get_seq_annotations(sequence)
        for cann in annotations:
            shortdesc.append(({'annotationtype': 'other', 'sequence': sequence}, cann))
        return shortdesc

    def show_annotation_info(self, annotation):
        '''Show the website for the sequence

        Parameters
        ----------
        annotation : dict
            should contain 'sequence'
        '''
        # open in a new tab, if possible
        new = 2

        address = '%s/search_results?sequence=%s' % (self.web_interface, annotation['sequence'])
        webbrowser.open(address, new=new)
