import pandas
from django.conf import settings

# setup django do be able to access django db models 
import GEN_database.settings as GEN_settings

try:
    settings.configure(INSTALLED_APPS=GEN_settings.INSTALLED_APPS,
                       DATABASES=GEN_settings.DATABASES)

    import django
    django.setup()
except:
    print("django setup failed-- already done?")
    pass

class DB:
    def __init__(self, db_path):
        import sqlite3
        
        self.db_path = db_path
        print("connecting to ", self.db_path)
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        
    def get_fastq_metadata(self, metric_name, index_str=True, analysis_id=False):
        
        if analysis_id:
            analysis_filter = f'and analysis_id={analysis_id}'
        else:
            analysis_filter = ''
        sql = f'''select fastq_id,value from  GEN_fastqfilesmetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  where t2.name="{metric_name}" {analysis_filter}            
                  union 
                  select t3.fastq_id,t1.value from  GEN_samplemetadata t1
                  inner join GEN_term t2 on t1.term_id=t2.id
                  inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
                  where t2.name="{metric_name}"
                  '''

        print(sql)
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_sample_metadata(self, metric_name, index_str=True):
        sql = f'''select t1.id,value from GEN_sample t1
                  inner join GEN_samplemetadata t2 on t1.id=t2.sample_id
                  inner join GEN_term t3 on t2.term_id=t3.id
                  where t3.name="{metric_name}";'''
                  
        if index_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_fastq_id2run_name(self,):

        sql = '''select t1.id,run_name from GEN_fastqfiles t1
                 inner join GEN_runs t2 on t1.run_id=t2.id
        '''

        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}

    def get_metadata_labels(self,):

        sql = '''select distinct t2.id,t2.name from GEN_fastqfilesmetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
                 union 
                 select distinct t2.id,t2.name from GEN_samplemetadata t1
                 inner join GEN_term t2 on t1.term_id=t2.id
        '''
        
        return {str(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_molis_id2fastq_id(self,):

        sql = ''' 
        select fastq_id, molis_id from GEN_fastqtosample t1
        inner join GEN_sample t2 on t1.sample_id=t2.id
        '''

        return {i[0]:i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_sample_df(self,):

        sql = ''' 
        select * from GEN_sample
        '''

        return pandas.read_sql(sql, self.conn)


    def get_fastq_metadata_list(self, 
                                term_list=False, 
                                fastq_filter=None, 
                                run_name_list=False,
                                metadata_value_list=False,):
        '''
        retrieve metadata from both sample and fastq metadata table
        '''

        res_filter = ''
        if term_list:
            term_filter = '","'.join(term_list)
            res_filter += f'and t2.name in ("{term_filter}")\n' 
        if fastq_filter:
            fastq_filter_str = ','.join([str(i) for i in fastq_filter])
            res_filter += f'and fastq_id in ({fastq_filter_str})\n' 
        if run_name_list:
            run_filter = '","'.join(run_name_list)
            res_filter += f'and run_name in ("{run_filter}")'
        if metadata_value_list:
            metadata_filter = '","'.join(metadata_value_list)
            res_filter += f'and t1.value in ("{metadata_filter}")'
        
        
        sql = f'''
            select distinct fastq_id, t2.name,t1.value, run_name from GEN_fastqfilesmetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqfiles t3 on t1.fastq_id=t3.id 
            inner join GEN_runs t4 on t3.run_id=t4.id 
            where t3.fastq_prefix not like "Undetermined%"
            {res_filter}
            group by fastq_id, t2.name,t3.fastq_prefix,t4.run_date,t4.run_name,t4.read_length
            union 
            select distinct t3.fastq_id, t2.name,t1.value,run_name from GEN_samplemetadata t1 
            inner join GEN_term t2 on t1.term_id=t2.id 
            inner join GEN_fastqtosample t3 on t1.sample_id=t3.sample_id
            inner join GEN_fastqfiles t4 on t3.fastq_id=t4.id 
            inner join GEN_runs t5 on t4.run_id=t5.id 
            where t4.fastq_prefix not like "Undetermined%"
            {res_filter}
            group by t3.fastq_id, t2.name,t4.fastq_prefix,t5.run_date,t5.run_name,t5.read_length
            '''
        print(sql)
        # AND t4.run_name like "%_CleanPlex"
        # 
        df = pandas.read_sql(sql, self.conn)


        return df

    def get_xslx_id2fastq_id(self,):
        
        sql = '''select t2.xlsx_sample_ID,t1.fastq_id from GEN_fastqtosample t1 
                inner join GEN_sample t2 on t1.sample_id=t2.id '''

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_xslx_id2sample_id(self,):
        
        sql = '''select xlsx_sample_ID,id from GEN_sample'''

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def count_qc_warning(self, key_str=False):
        
        df = self.get_fastq_metadata_list()[["fastq_id","value"]]

        term2n_warnings = df.query('value == "WARN"').groupby(["fastq_id"])['value'].count().to_dict()

        if key_str:
            return {str(i):term2n_warnings[i] for i in term2n_warnings}
        else:
            return term2n_warnings




    def get_fastq_metadata_stats(self, term_list):

        df = self.get_fastq_metadata_list(term_list)[["fastq_id","name","value"]]
        
        df["value"] = pandas.to_numeric(df["value"])
        term2median = df.groupby(["name"])['value'].median().round(3).to_dict()
        term2mean = df.groupby(["name"])['value'].mean().round(3).to_dict()
        term2sd = df.groupby(["name"])['value'].std().round(3).to_dict()        

        return term2median, term2mean, term2sd

    def get_fastq_and_sample_data(self, fastq_id_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_id_list])
        
        # left join because some fastq won't have match in the sample table
        # possible problem: fastq prefix match with multiple samples from different species
        # in that case: remove species name
        sql = f'''select distinct t1.id as fastq_id,fastq_prefix,R1,R2,species_name,molis_id from GEN_fastqfiles t1 
                left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
                left join GEN_sample t3 on t2.sample_id=t3.id 
                where t1.id in ("{fastq_list_filter}");
            '''
        print(sql,)
        return pandas.read_sql(sql, self.conn)


    def get_analysis_metadata(self, term_name):
        
        sql = f"""select t1.analysis_id,value from GEN_analysismetadata t1 
                  inner join GEN_term t2 on t1.term_id =t2.id 
                  where t2.name like '{term_name}';"""

        return {int(i[0]):i[1] for i in self.cursor.execute(sql,).fetchall()}


    def get_run_name2run_id(self,):
        sql = 'select run_name,id from GEN_runs'
        return {i[0]:i[1] for i in self.cursor.execute(sql,).fetchall()}

    def match_fastq_to_sample(self, fastq_prefix):

        # mapping based on xlsx sample id
        # conflict possible with several projects with sample with numeric labels (1,2,3, 232,...)
        # use date filter to exclude mapping older than X on the various columns
        sql = f'select id from GEN_sample where xlsx_sample_ID="{fastq_prefix}"' 
        try:
            sample_id = self.cursor.execute(sql,).fetchall()[0][0]
        except:
            sql = f'select id from GEN_sample where sample_name="{fastq_prefix}"' 
            try:
                sample_id = self.cursor.execute(sql,).fetchall()[0][0]
            except:
                sql = f'select id from GEN_sample where alias="{fastq_prefix}"' 
                try:
                    sample_id = self.cursor.execute(sql,).fetchall()[0][0]
                except:
                    sample_id = None 
        return sample_id

    def get_sample_id(self, sample_xls_id):
        sql = 'select id from GEN_sample where xlsx_sample_ID=?'

        return self.cursor.execute(sql,(sample_xls_id,)).fetchall()[0][0]

    def add_sample_to_fastq_relation(self, fastq_id, sample_id):
        
            # ignore if relation already known
            sql2 = 'insert or ignore into GEN_fastqtosample(fastq_id, sample_id) values(?,?)'
            
            self.cursor.execute(sql2, 
                               (fastq_id,sample_id))
            
            self.conn.commit()


    def insert_run(self,
                   run_name,
                   run_date,
                   assay,
                   read_length,
                   paired,
                   filearc_folder):
        
        sql = '''INSERT into GEN_runs (run_name, run_date, assay, read_length, paired, filearc_folder) values(?,?,?,?,?,?) 
            ON CONFLICT(GEN_runs.run_name) DO UPDATE SET run_date=?, assay=?, read_length=?, paired=?, filearc_folder=?;
        '''
        self.cursor.execute(sql, [run_name,
                                  run_date,
                                  assay,
                                  read_length,
                                  paired, 
                                  filearc_folder,
                                  run_date,
                                  assay,
                                  read_length,
                                  paired, 
                                  filearc_folder])
        self.conn.commit()

    def insert_sample(self,
                      col_names,
                      values_list,
                      sample_xls_id):
        
        update_str = '%s=?'

        update_str_comb = ','.join([update_str % colname for colname in col_names])
        # INSERT into GEN_sample(xlsx_sample_ID,species_name,date_received,sample_name,sample_type,analysis_type,description,molis_id,myseq_passage,run_date,date_registered,date_sample_modification,user_creation_id,user_modification_id) values(?,?,?,?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET xlsx_sample_ID=?,species_name=?,date_received=?,sample_name=?,sample_type=?,analysis_type=?,description=?,molis_id=?,myseq_passage=?,run_date=?,date_registered=?,date_sample_modification=?,user_creation_id=?,user_modification_id=?;
        # [29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2, 29, 'Staphylococcus aureus', False, nan, 'strain', 'research', nan, 1306182530, nan, False, '2020-09-02', '2020-09-02', 2, 2]            
        # update all columns in case of conflict with xlsx_sample_id
        
        # NOTE: xlsx_sample_ID used as reference: if a row is updated in the xlsx table, the corresponding row is updated in the sql table
        sql_template = 'INSERT into GEN_sample(%s) values(%s)' \
                       ' ON CONFLICT(GEN_sample.xlsx_sample_ID) DO UPDATE SET %s;' % (','.join(col_names),
                                                                                      ','.join(['?']*len(col_names)),
                                                                                      update_str_comb)
        
        print(sql_template)
        print(values_list)
        self.cursor.execute(sql_template, values_list + values_list)
        self.conn.commit()

        return self.get_sample_id(sample_xls_id)
  

    def match_sample_to_fastq(self, sample_prefix):
        sql = 'select id from GEN_fastqfiles where fastq_prefix=?' 
        try:
            fastq_id_list = [i[0] for i in self.cursor.execute(sql,(sample_prefix,)).fetchall()]
        except:
            fastq_id = [] 
        return fastq_id_list


    def get_fastq(self, run_name=False):
        
        sql = 'select t1.id,run_name,date_run,qc,fastq_prefix,xlsx_sample_ID,species_name,date_received,read_length,t1.id from GEN_fastqfiles t1 ' \
            ' inner join GEN_runs t2 on t1.run_id=t2.id ' \
            ' left join GEN_fastqtosample t3 on t1.id=t3.fastq_id' \
            ' left join GEN_sample t4 on t3.sample_id=t4.id'

        sql = '''
        select t1.id as fastq_id, run_name,t5.status,fastq_prefix,read_length,t3.id as sample_id,taxonomy,date_received,t4.run_date,t1.id,t4.qc_id from GEN_fastqfiles t1 
        left join GEN_fastqtosample t2 on t1.id=t2.fastq_id 
        left join GEN_sample t3 on t2.sample_id=t3.id 
        left join GEN_runs t4 on t1.run_id=t4.id
        left join GEN_analysis t5 on t4.qc_id=t5.id
        '''

        if run_name:
            sql += f'\nwhere run_name="{run_name}"'

        #print(sql)
    
        return pandas.read_sql(sql, self.conn)
      
    
    def get_run_table(self,):
        
        from GEN.models import Analysis
        
        sql = '''select run_date,run_name,read_length,filearc_folder,qc_id,qc_path,count(*) as n_fastq from GEN_runs t1
        inner join GEN_fastqfiles t2 on t1.id=t2.run_id group by run_date,run_name,read_length,filearc_folder,qc_id,qc_path
        ''' 
        data = [list(i) for i in self.cursor.execute(sql,).fetchall()]
        for n, row in enumerate(data):
            if row[4]:
                data[n][4] = Analysis.objects.filter(id=row[4])[0]
        return data

    def get_run_samples(self, run_name):
        
        sql = f'''select fastq_prefix,R1,R2 from GEN_fastqfiles t1 
                  inner join GEN_runs t2 on t1.run_id=t2.id where run_name="{run_name}"
               '''
        return self.cursor.execute(sql,).fetchall()
    
    
    def get_run_sample2species(self, run_name):
        # left join because some fastq won't have match in the sample table
        sql = f'''select fastq_prefix,species_name from GEN_fastqfiles t1 
                inner join GEN_runs t2 on t1.run_id=t2.id
                left join GEN_fastqtosample t3 on t1.id=t3.fastq_id
                left join GEN_sample t4 on t3.sample_id=t4.id where run_name="{run_name}";
            '''
        return {i[0]: i[1] for i in self.cursor.execute(sql,).fetchall()}
    
        
    def fastq_prefix2fastq_id(self,fastq_prefix_list):
        
        fasta_filter = '","'.join(fastq_prefix_list)
        sql = f'select fastq_prefix,id from GEN_fastqfiles where fastq_prefix in ("{fasta_filter}")'
        
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def add_fastq_metadata(self, 
                           fastq_id,
                           term_name,
                           value,
                           analysis_id=False):
        '''
        [{"fastq_id": fastq_id,
        "term_name": <name>,
        "value": <value>,
        "analysis_id" analysis_id}]
        '''
        from GEN.models import Term
        from GEN.models import FastqFilesMetadata
        term = Term.objects.get_or_create(name=term_name)[0]
        if not analysis_id:
            m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value)
            m.save()
        else:
            m = FastqFilesMetadata(term=term, fastq_id=fastq_id, value=value, analysis_id=analysis_id)
            m.save()

    def add_sample_metadata(self, 
                           sample_id,
                           term_name,
                           value):
        '''
        [{"sample_id": sample_id,
        "term_name": <name>,
        "value": <value>,
        "analysis_id" analysis_id}]
        '''
        from GEN.models import Term
        from GEN.models import SampleMetadata
        term = Term.objects.get_or_create(name=term_name)[0]

        m = SampleMetadata(term=term, sample_id=sample_id, value=value)
        m.save()


    def get_term2term_id(self, 
                         term_list):
        from GEN.models import Term
        term2term_id = {}
        for term_name in term_list:
            term = Term.objects.get_or_create(name=term_name)[0]
            term2term_id[term_name] = term.id
        return term2term_id

    def get_term_id2term_name(self, 
                              term_id_list):
        from GEN.models import Term
        term_id2term_name = {}
        for term_id in term_id_list:
            term = Term.objects.get_or_create(id=term_id)[0]
            term_id2term_name[term_id] = term.name
        return term_id2term_name

    def fastq_mutation(self, aa_change):
        
        sql = f'select fastq_id from GEN_snps where aa_change="{aa_change}";'

        return [i[0] for i in self.cursor.execute(sql,)]


    def add_QC_report(self, run_name, run_path):
        
        sql = 'update GEN_runs set qc=1 where run_name=?'
        self.cursor.execute(sql, [run_name]) 
        
        sql = 'update GEN_runs set qc_path=? where run_name=?'
        self.cursor.execute(sql, (run_path, run_name))
        self.conn.commit()
    
    
    def fastq_qc_filter(self, analysis_id, value):
        
        sql = f'''select fastq_id from GEN_fastqfilesmetadata t1 
                 inner join GEN_term t2 on t1.term_id=t2.id 
                 where t1.analysis_id={analysis_id} and t2.name="qc_status" and t1.value="{value}";
        '''
        print(sql)
        
        return [i[0] for i in self.cursor.execute(sql,)]
        
    def format_snps(self, fastq_id_list):
        
        fastq_filter = ','.join([str(i) for i in fastq_id_list])
        sql = f'select fastq_id,nucl_change, aa_change, gene from GEN_snps where fastq_id in ({fastq_filter})'

        df = pandas.read_sql(sql, self.conn).set_index("fastq_id")
        
        print("head", df.head())
        
        fastq2data = {fastq_id: {} for fastq_id in df.index.unique()}
        
        for fastq_id in list(fastq2data.keys()):
            target_fastq = df.loc[fastq_id]
            aa_changes = [f'{pos["aa_change"]} ({pos["gene"]})' for n, pos in target_fastq.iterrows() if not pandas.isna(pos["aa_change"])]
            nucl_changes = [pos["nucl_change"] for n, pos in target_fastq.iterrows()]
            fastq2data[fastq_id]["aa_changes"] = ';'.join(aa_changes)
            fastq2data[fastq_id]["nucl_changes"] = ';'.join(nucl_changes)
            
        return fastq2data   
        
        
        
    
    def insert_or_get_fastq(self, 
                            fastq_prefix, 
                            run_id, 
                            R1, 
                            R2):
        from GEN.models import FastqFiles

        if "control" in fastq_prefix:
            control = 1
        else:
            control = 0
        # insert 
        fastq = FastqFiles.objects.get_or_create(control_sample=control,
                                                 fastq_prefix=fastq_prefix,
                                                 run_id=run_id,
                                                 R1=R1,
                                                 R2=R2)[0]
        
        return fastq.id

    def get_sample2species(self, sample_list):
        
        sample_list_filter = '","'.join(sample_list)
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where fastq_prefix in ("{sample_list_filter}");
           '''
        
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def get_fastq_id2species(self, fastq_list):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_prefix,species_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        print(sql)
        return {i[0]:i[1] for i in self.cursor.execute(sql,)}
    
    def get_fastq_id2sample_name(self, fastq_list, key_str=True):
        
        fastq_list_filter = '","'.join([str(i) for i in fastq_list])
        sql = f'''select distinct fastq_id,sample_name from GEN_fastqfiles t1 
              left join GEN_fastqtosample t2 on t1.id=t2.fastq_id
              left join GEN_sample t3 on t2.sample_id=t3.id 
              where t1.id in ("{fastq_list_filter}");
           '''
        print(sql)
        if key_str:
            return {str(i[0]):i[1] for i in self.cursor.execute(sql,)}
        else:
            return {int(i[0]):i[1] for i in self.cursor.execute(sql,)}

    def get_sample_table(self,):
        
        sql = '''
        select distinct A.id,A.sample_type,A.sample_name,A.taxonomy,A.date_registered, A.date_received,A.n_fastq, count(t3.subproject_id ) as n_projects from (
        select distinct t1.id,t1.sample_type,t1.sample_name,t1.taxonomy,t1.date_registered, t1.date_received, count(t2.fastq_id) as n_fastq from GEN_sample t1 
        left join GEN_fastqtosample t2 on t1.id=t2.sample_id 
        group by t1.id) A  
        left join GEN_subprojectsample t3 on A.id=t3.sample_id
        group by A.id 
        '''

        return self.cursor.execute(sql,).fetchall()
    
    def get_analysis_fastq_list(self,analisis_id):
        
        sql = f'''
        select fastq_id from GEN_fastqset where analysis_id={analisis_id}
        '''
        
        return [i[0] for i in self.cursor.execute(sql,).fetchall()]
        

def update_analysis_status(analysis_id, status):
    from GEN.models import Analysis
    m = Analysis.objects.filter(id=analysis_id)[0]
    m.status = status
    m.save()
    
def add_analysis_metadata(analysis_id, term, value, update=False):
    from GEN.models import Term
    from GEN.models import AnalysisMetadata
    term = Term.objects.get_or_create(name=term)[0]
    if not update:
        m = AnalysisMetadata(term=term, analysis_id=analysis_id, value=value)
        m.save()
    else:
        # update value
        m = AnalysisMetadata.objects.filter(term=term, analysis_id=analysis_id)[0]
        m.value = value
        m.save()
        
def create_analysis(fastq_id_list,
                    analysis_description,
                    reference_fastq_id_list = [],
                    subproject_id=False,
                    workflow_name="Airflow_epidemiology"):
                
    from GEN.models import Workflow
    from GEN.models import WorkflowSteps
    from GEN.models import FastqFiles
    from GEN.models import FastqSet
    from GEN.models import Term
    from GEN.models import AnalysisStatus
    from GEN.models import ProjectAnalysis
    from GEN.models import Analysis
    from GEN.models import AnalysisMetadata
    from datetime import datetime
    
    # create new analysis
    workflow_instance = Workflow.objects.filter(workflow_name=workflow_name)
    workflow_instance = workflow_instance[0]
    analysis = Analysis(workflow=workflow_instance, start_date=datetime.today(), status='started')
    analysis.save()
    
    # associated to project if a project id was provided
    if subproject_id:
        project_analysis = ProjectAnalysis(analysis=analysis, subproject_id=subproject_id)
        project_analysis.save()
    
    # insert description as metadata
    term = Term.objects.get_or_create(name="description")[0]
    desc = AnalysisMetadata(term=term, analysis=analysis, value=analysis_description)
    desc.save()
    
    if len(reference_fastq_id_list) > 0:
        term_ref_genome = Term.objects.get_or_create(name="reference_genome")[0]
        for ref_genome in reference_fastq_id_list:
            m = AnalysisMetadata(term=term_ref_genome, analysis=analysis, value=ref_genome)
            m.save()  
    
    # add each fastq to fastq set
    for fastq_id in fastq_id_list:
        print("fastq_id", fastq_id)
        fastq_instance = FastqFiles.objects.filter(id=fastq_id)[0]
        new_set_fastq = FastqSet(analysis=analysis, fastq=fastq_instance)
        new_set_fastq.save()

    # create "AnalysisStatus" entry for each step of the workflow
    # ==> all marked as not DONE
    all_steps = WorkflowSteps.objects.filter(workflow=workflow_instance)
    for one_step in all_steps:
        new_project_analysis_status_instance = AnalysisStatus(analysis=analysis, step=one_step.step, status=0)
        new_project_analysis_status_instance.save()
        
    
    return analysis.id
