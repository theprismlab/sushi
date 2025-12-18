import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

String sectionHeaderStyleGreen = ' color: white; background: #dbdb8e; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 30px '
String sectionHeaderStyleBlue = ' color: white; background: #7ea6d3; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 18px'
String sectionHeaderStyleRed = ' color: white; background: red; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 30px'
String sectionHeaderStyleYellow = ' color: white; background: yellow; font-family: Roboto, sans-serif !important; padding: 5px; text-align: center; font-size: 12px'
String separatorStyleCss = ' border: 0; border-bottom: 1px dashed #ccc; background: #999; '

pipeline {
    agent any
    // Define parameters that can be edited via the Jenkins UI
    parameters {
        separator(
          name: "user_inputs",
          sectionHeader: "User Inputs",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleGreen
        )

        separator(
          name: "run_sushi",
          sectionHeader: "Run Sushi",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'TRIGGER_BUILD', defaultValue: true, description: 'Check this to trigger the build jobs below. If unchecked, the build will not be triggered and only the config.json will be generated.')

        separator(
          name: "upload_to_bq",
          sectionHeader: "Upload results to BigQuery",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'UPLOAD_TO_BQ', defaultValue: false, description: 'Check this to upload the results to BigQuery.')

        separator(
          name: "build_details",
          sectionHeader: "Build Details",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        string(name: 'BUILD_DIR', defaultValue: '/cmap/obelix/pod/prismSeq/', description: 'Directory where the build output will go. Must contain the raw counts file from Nori. If not getting your metadata from cellDB & COMET this directory must also include the sample and cell line/pool metadata.')
        string(name: 'BUILD_NAME', defaultValue: '', description: 'Build name; used to name output files from the adapter and QC scripts. Should match the directoty name from BUILD_DIR.')
        string(name: 'SCREEN', defaultValue: '', description: 'Name of the screen with which this build is associated. This is necessary if you are getting your metadata from COMET/cellDB and/or you plan to upload the data to the portal. This should be the same as the screen from the sample metadata file.')
        choice(name: 'SCREEN_TYPE', choices: ['MTS_SEQ', 'CPS_SEQ', 'EPS_SEQ', 'APS_SEQ'], description: 'Select the type of build you are running.')
        string(name: 'PERT_PLATES', defaultValue: '', description: 'Comma separated list of pert_plates to include in this build. If not provided, all plates for the given screen will be used.')
        string(name: 'SIG_COLS', defaultValue: 'cell_set,pert_type,pert_name,pert_id,pert_dose,pert_dose_unit,day,x_project_id,pert_plate,pert_vehicle', description: 'List of signature columns found in the sample meta that describe unique treatment conditions. Generally, this list should NOT include replicate information such as \"tech_rep\" or \"bio_rep\".')

        separator(
          name: "metadata",
          sectionHeader: "Metadata",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'CREATE_CELLDB_METADATA', defaultValue: true, description: 'Get cell line and pool metadata from cellDB. Requires that all sample_meta column entries for cell_set exist in cellDB.')
        booleanParam(name: 'CREATE_SAMPLE_META', defaultValue: false, description: 'Get sample metadata from COMET, project must be registered and all metadata steps completed.')

        separator(
            name: "controls",
            sectionHeader: "Controls",
            separatorStyle: separatorStyleCss,
            sectionHeaderStyle: sectionHeaderStyleBlue
        )
        string(name: 'CONTROL_BARCODE_META', defaultValue: 'h-a', description: 'Metadata for control barcodes. If the CBs exist in cellDB, this can simply be the lowercase cb_ladder name (ie, h-a). Otherwise, this must be a csv file located in the build directory.')
        choice(name: 'NORM_METHOD', choices: ['normalize_then_shift', 'normalize_with_pseudocount'], description: 'Select the normalization method to use.')
        string(name: 'CTL_TYPES', defaultValue: 'ctl_vehicle', description: 'Value in the pert_type column of the sample meta that identifies the negative contols.')
        string(name: 'POSCON_TYPE', defaultValue: 'trt_poscon', description: 'Value in the pert_type column of the sample meta that identifies the positive controls.')
        string(name: 'CONTROL_COLS', defaultValue: 'cell_set,day,pcr_plate,replicate_plate,pert_vehicle', description: 'List of columns found in the sample meta that describe individual negative control conditions.')

        separator(
          name: "core_modules",
          sectionHeader: "Core Modules",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'COLLATE_FASTQ_READS', defaultValue: true, description: 'Checks to ensure raw reads come from expected flowcells and lanes and then sums the counts across samples (SEQUENCING_INDEX_COLS).')
        booleanParam(name: 'FILTER_COUNTS', defaultValue: true, description: 'Assigns raw reads to the appropriate treatment conditions and cell lines; filters those that do not match. Removes cell lines that are duplicated in a cell set.')
        booleanParam(name: 'FILTER_SKIPPED_WELLS', defaultValue: true, description: 'Check this to filter out wells that were skipped by the echo.')
        booleanParam(name: 'REMOVE_DATA', defaultValue: false, description: 'Uses a data_to_remove.csv files to remove data. Runs as part of filter counts.')
        booleanParam(name: 'CBNORMALIZE', defaultValue: true, description: 'Normalizes counts. Requires vehicle controls and a control barcode ladder.')
        booleanParam(name: 'COMPUTE_LFC', defaultValue: true, description: 'Compute the fold changes from vehicle controls of each cell line for each treatment condition.')
        booleanParam(name: 'BIAS_CORRECTION', defaultValue: true, description: 'Correct fold change values using cell line growth annotations.')
        booleanParam(name: 'COLLAPSE', defaultValue: true, description: 'Median collapses biological replicates.')

        separator(
          name: "qc_modules",
          sectionHeader: "QC Modules",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'GENERATE_QC_TABLES', defaultValue: true, description: 'Generate MTS style QC tables')
        booleanParam(name: 'GENERATE_QC_TABLES_2', defaultValue: false, description: 'Generate QC tables post-lfc. Only check this if you are computing log fold changes.')
        booleanParam(name: 'QC_IMAGES', defaultValue: false, description: 'Check this to trigger the QC images job.')
        booleanParam(name: 'FILTER_QC_FLAGS', defaultValue: true, description: 'Filters data that have qc_flags.')
        booleanParam(name: 'FILTER_FAILED_LINES', defaultValue: true, description: 'Filter cell lines from LFC and collapsed LFC that have failed QC.')

        separator(
          name: "portal_prep",
          sectionHeader: "Portal Prep",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'CONVERT_SUSHI', defaultValue: false, description: 'Convert output column headers to format for MTS pipeline and upload to s3.')

        separator(
          name: "analytics_modules",
          sectionHeader: "Analytics Modules",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleBlue
        )
        booleanParam(name: 'DRC', defaultValue: false, description: 'Generate dose response curves.')
        booleanParam(name: 'UNIVARIATE_BIOMARKER', defaultValue: false, description: 'Run univariate biomarker analysis.')
        booleanParam(name: 'MULTIVARIATE_BIOMARKER', defaultValue: false, description: 'Run multivariate biomarker analysis.')
        booleanParam(name: 'LFC_BIOMARKER', defaultValue: false, description: 'Use log fold change values to run biomarker analysis. Requires the COMPUTE_LFC module to have been run.')
        booleanParam(name: 'AUC_BIOMARKER', defaultValue: false, description: 'Use AUC values to run biomarker analysis. Requires the DRC module to have been run.')
        booleanParam(name: 'SYNERGY', defaultValue: false, description: 'Compute synergy for CPS.')
        booleanParam(name: 'SYNERGY_BIOMARKER', defaultValue: false, description: 'Use synergy values to run biomarker analysis. SYNERGY module to run.')

        // Parameters that we don't expect users to change
        separator(
          name: "do_not_edit",
          sectionHeader: "WARNING: ADVANCED USERS ONLY. DO NOT EDIT",
          separatorStyle: separatorStyleCss,
          sectionHeaderStyle: sectionHeaderStyleRed
        )

        // pipeline version
        string(name: 'GIT_BRANCH', defaultValue: 'main', description: 'Pipeline branch to use')
        booleanParam(name: 'USE_LATEST', defaultValue: true, description: 'Check this to use the most up to date version from the specified branch. If not checked, will use the specified commit.')
        string(name: 'COMMIT_ID', defaultValue: '', description: 'Specific commit ID to use (leave empty if using the latest commit in the branch or if already specified in the config file.)')

        // QC Paramters
        string(name: 'QC_PARAMS', defaultValue: 'qc_params.json', description: 'File name in BUILD_DIR containing the QC parameters.')
        string(name: 'nc_variability_threshold', defaultValue: '1', description: 'Threshold for negative control variability')
        string(name: 'error_rate_threshold', defaultValue: '0.05', description: 'Threshold for error rate')
        string(name: 'pc_viability_threshold', defaultValue: '0.25', description: 'Threshold for positive control viability')
        string(name: 'nc_raw_count_threshold', defaultValue: '40', description: 'Threshold for negative control raw counts')
        string(name: 'contamination_threshold', defaultValue: '0.8', description: 'Threshold for fraction of expected reads')
        string(name: 'cb_mae_threshold', defaultValue: '1', description: 'Threshold for mean absolute error of control barcodes')
        string(name: 'cb_spearman_threshold', defaultValue: '0.8', description: 'Threshold for control barcode spearman correlation')
        string(name: 'cb_cl_ratio_low_negcon', defaultValue: '0', description: 'Low threshold for control barcode ratio in negative controls')
        string(name: 'cb_cl_ratio_high_negcon', defaultValue: '100', description: 'High threshold for control barcode ratio in negative controls')
        string(name: 'cb_cl_ratio_low_poscon', defaultValue: '0', description: 'Low threshold for control barcode ratio in positive controls')
        string(name: 'cb_cl_ratio_high_poscon', defaultValue: '100', description: 'High threshold for control barcode ratio in positive controls')
        string(name: 'well_reads_threshold', defaultValue: '40', description: 'Minimum median control barcode reads per well')
        string(name: 'pool_well_delta_threshold', defaultValue: '5', description: 'Maximum delta of log2_normalized_n between a cell line and the pool median in a given well before it is considered an outlier')
        string(name: 'pool_well_fraction_threshold', defaultValue: '0.4', description: 'Minimum fraction of cells in a pool that must be outliers in order to flag that pool/well')
        string(name: 'fraction_expected_controls', defaultValue: '0.667', description: 'Fraction of expected controls that must be present in a given pcr_plate for it to be considered a valid well. If either vehicle or poscon wells fall below this threshold, the entire pcr_plate will be removed.')
        // Sequencing tech
        choice(name: 'SEQ_TYPE', choices: ['DRAGEN', 'MiSeq', 'HiSeq', 'NovaSeq'], description: 'Choose DRAGEN, MiSeq, HiSeq, or NovaSeq. MiSeq and HiSeq/NovaSeq return files named differently. This setting sets the INDEX_1, INDEX_2, and BARCODE_SUFFIX parameters in fastq2readcount. Select DRAGEN if fastq files are from the DRAGEN pipeline from GP. Choosing NovaSeq reverses index 2.')

        // Sushi Input files
        string(name: 'RAW_COUNTS_UNCOLLAPSED', defaultValue: 'raw_counts_uncollapsed.csv.gz', description: 'File name in BUILD_DIR containing the uncollapsed raw counts. This should be the file generated by Nori.')
        string(name: 'SAMPLE_META', defaultValue: 'sample_meta.csv', description: 'File name in BUILD_DIR of the sample meta.')
        string(name: 'CELL_SET_AND_POOL_META', defaultValue: 'cell_set_and_pool_meta.csv', description: 'Cell set and pool information for this run.')
        string(name: 'CELL_LINE_META', defaultValue: 'cell_line_meta.csv', description: 'File in BUILD_DIR containing cell line metadata')
       
        // Additional parameters ordered by when they first appear
        // Collate FASTQ reads
        string(name: 'BARCODE_COL', defaultValue: 'forward_read_barcode', description: 'Name of the column containing the barcode sequence. The column containing the barcode sequence should have the same name across the Nori output file, the cell line metadata, and the CB metadata. This defaults to \"forward_read_barcode\", and the paramter is first used in COLLATE_FASTQ_READS.')
        string(name: 'SEQUENCING_INDEX_COLS', defaultValue: 'flowcell_names,index_1,index_2', description: 'List of sequencing related columns found in the sample meta. These columns are used to uniquely identify every PCR well in the run. Default value is \"flowcell_names,index_1,index_2\", and the parameter is used in COLLATE_FASTQ_READS.')
        string(name: 'ID_COLS', defaultValue: 'pcr_plate,pcr_well', description: 'List of columns found in the sample meta that are used to create a unique ID for each sample-replicate. This defaults to \"pcr_plate,pcr_well\", but could be any combination sample meta columns that uniquely identifies every well in the run. This parameter is first used in COLLATE_FASTQ_READS.')
        string(name: 'CHUNK_SIZE', defaultValue: '10000000', description: 'Number of rows for a chunk. Due to the large size of the Nori output, some actions are performed in chunks to conserve memory. This parameter sets the size of a chunk and defaults to 10^6 or \"10000000\". This paramter is first used in COLLATE_FASTQ_READS.')
        string(name: 'LOW_ABUNDANCE_THRESHOLD', defaultValue: '20', description: 'Threshold for unknown barcodes. Unknown barcodes below this threshold will be deidentified and their counts will be included under the term unknown_low_abundance_barcode. This paramter defaults to \"20\" and is used in COLLATE_FASTQ_READS.')
        // Normalize
        string(name: 'PSEUDOCOUNT', defaultValue: '0', description: 'Pseudocount value added to all reads before log transformations. This defaults to \"0\" and is used in CBNORMALIZE.')
        string(name: 'READ_DETECTION_LIMIT', defaultValue: '10', description: 'Smallest read count value used to compute pseudovalues.')
        // Compute l2fc 
        string(name: 'CELL_LINE_COLS', defaultValue: 'pool_id,depmap_id,lua,cell_set,growth_condition', description: 'List of columns across the metadata files that are used to identify a unique cell line. This defaults to \"pool_id,depmap_id,lua\", but can also include \"cell_set\" or descriptive columns like \"project_code\" that you would like to pass through the pipeline. This parameter is first used in COMPUTE_LFC.')
        string(name: 'COUNT_COL_NAME', defaultValue: 'log2_normalized_n', description: 'Name of the numerical column that should be used to compute log2 fold change values. This defaults to \"normalized_n\" and is used in COMPUTE_LFC.')
        string(name: 'COUNT_THRESHOLD', defaultValue: '40', description: 'Threshold for filtering the negative controls. In the negative control conditions, cell lines whose median counts are below this threshold are not confidently detected and thus are dropped. This defaults to \"40\" and is used in COMPUTE_LFC.')
        string(name: 'BIO_REP_COL', defaultValue: 'bio_rep', description: 'Column identifying the biological replicates. Defaults to \"bio_rep\".')
        // Collapse replicates
        string(name: 'L2FC_COLUMN', defaultValue: 'l2fc', description: 'Name of the column containing the log2 fold change values used in DRC. This defaults to \"l2fc\".')
        string(name: 'GROWTH_PATTERN_COL', defaultValue: 'growth_condition', description: 'Name of the column containing the cell line growth annotations. This defaults to \"growth_condition\".')
        string(name: 'COLLAPSED_L2FC_COLUMN', defaultValue: 'median_l2fc', description: 'Name of the column containing the collapsed log2 fold change values used in biomarker. This defaults to \"collapsed_l2fc\".')
        // DRC
        string(name: 'VIABILITY_CAP', defaultValue: '1.5', description: 'Cap for viability values used when computing LFC. This defaults to \"1.5\".')
        // CPS
        string(name: 'COMBINATION_COL', defaultValue: 'is_combination', description: 'Name of the column describing if a row has a combination. This defaults to \"is_combination\".')
        string(name: 'N_SAMPLES', defaultValue: '10000', description: 'Size of the reampling. This defaults to \"10000\".')
        string(name: 'SYNERGY_L2FC_COL', defaultValue: 'median_l2fc_uncorrected', description: '[SYNERGY] Column name in collapsed_l2fc.csv to be used for synergy calculations. This defaults to \"median_l2fc_uncorrected\".')

        // Files created by sushi
        string(name: 'PRISM_BARCODE_COUNTS', defaultValue: 'prism_barcode_counts.csv', description: 'Filename in BUILD_DIR containing PRISM barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'UNKNOWN_BARCODE_COUNTS', defaultValue: 'unknown_barcode_counts.csv', description: 'Filename in BUILD_DIR containing unknown barcode counts. This file is created by COLLATE_FASTQ_READS.')
        string(name: 'ANNOTATED_COUNTS', defaultValue: 'annotated_counts.csv', description: 'Filename in BUILD_DIR containing annotated counts. This file is creaed by FILTER_COUNTS.')
        string(name: 'FILTERED_COUNTS', defaultValue: 'filtered_counts.csv', description: 'Filename in BUILD_DIR containing filtered counts. This file is created by FILTER_COUNTS.')
        string(name: 'NORMALIZED_COUNTS', defaultValue: 'normalized_counts.csv', description: 'Filename in BUILD_DIR containing normalized counts. This file is created by CBNORMALIZE.')
        string(name: 'LFC', defaultValue: 'l2fc.csv', description: 'Filename containing log2 fold change values. This file is created by COMPUTE_LFC.')
        string(name: 'COLLAPSED_LFC', defaultValue: 'collapsed_l2fc.csv', description: 'Filename in BUILD_DIR containing replicate collapsed l2fc values. This file is created by COLLAPSED_LFC.')
        string(name: 'SKIPPED_WELLS', defaultValue: 'skipped_wells.csv', description: 'Filename in BUILD_DIR containing skipped wells. This file is created by create_sample_meta.')
        // Other
        string(name: 'API_URL', defaultValue: 'https://api.clue.io/api/', description: 'API URL')
        string(name: 'MERGE_PATTERNS', defaultValue: 'log2_auc_multivariate_biomarkers*,log2_auc_univariate_biomarkers*,median_l2fc_multivariate_biomarkers*,median_l2fc_univariate_biomarkers*,DRC_TABLE*', description: 'Patterns to search for when merging files by project. May be changed based on modules run.')

        // Biomarker
        string(name: 'BIOMARKER_FILE', defaultValue: '/cmap/obelix/pod/biomarker/depmap_public_24q4.h5', description: 'Biomarker reference file.')
        string(name: 'DR_COLUMN', defaultValue: 'log2_auc', description: 'Name of the column containing AUC values used in biomarker analysis.')
        string(name: 'DR_PATH', defaultValue: 'DRC_TABLE.csv', description: 'File in drc/BUILD_DIR containing dose response curve data. This file is created by DRC.')

        //Synergy for combinations
        string(name: 'SYNERGY_FILE', defaultValue: 'synergy_scores.csv', description: 'File in BUILD_DIR containing the synergy scores.')
        string(name: 'SYNERGY_COL', defaultValue: 'synergy', 
               description: 'Name of the column in SYNERGY_FILE containing the synergy scores.')
    }

    environment {
        CONFIG_FILE_PATH = "${params.BUILD_DIR}/config.json"
        QC_PARAMS_FILE_PATH = "${params.BUILD_DIR}/qc_params.json"
        GOOGLE_APPLICATION_CREDENTIALS = '/local/jenkins/.gcp_credentials.json'
    }

    stages {
        stage('Checkout') {
            steps {
                script {
                    def config = [:]
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    if (params.USE_LATEST) {
                        // Checkout the latest commit from the specified branch
                        checkout([$class: 'GitSCM',
                                  branches: [[name: "*/${params.GIT_BRANCH}"]],
                                  doGenerateSubmoduleConfigurations: false,
                                  extensions: [],
                                  userRemoteConfigs: scm.userRemoteConfigs
                        ])
                        // Overwrite the commit ID in the config with the latest commit
                        def latestCommitID = sh(script: 'git rev-parse --short HEAD', returnStdout: true).trim()
                        config.COMMIT = latestCommitID
                        echo "Using latest commit: ${latestCommitID}"
                    } else {
                        // Use the commit ID specified in the config.json
                        if (config.COMMIT) {
                            checkout([$class: 'GitSCM',
                                      branches: [[name: config.COMMIT]],
                                      doGenerateSubmoduleConfigurations: false,
                                      extensions: [[$class: 'RelativeTargetDirectory', relativeTargetDir: '']],
                                      submoduleCfg: [],
                                      userRemoteConfigs: scm.userRemoteConfigs
                            ])
                            echo "Using commit from config.json: ${config.COMMIT}"
                        } else {
                            error("COMMIT not specified in config.json and USE_LATEST is false.")
                        }
                    }
                }
            }
        }

        stage('Generate JSON Config') {
            steps {
                script {
                    def paramList = [
                        'SEQ_TYPE', 'API_URL', 'BUILD_DIR', 'INDEX_1', 'INDEX_2', 'BARCODE_SUFFIX', 'CREATE_CELLDB_METADATA',
                        'BUILD_NAME', 'CONVERT_SUSHI', 'REMOVE_DATA', 'FILTER_SKIPPED_WELLS',
                        'COUNTS', 'SCREEN', 'GENERATE_QC_TABLES', 'POSCON_TYPE', 'DRC', 'L2FC_COLUMN','COLLAPSED_L2FC_COLUMN',
                        'SKIPPED_WELLS','FILTER_QC_FLAGS', 'PERT_PLATES', 'SCREEN_TYPE',

                        // sushi input files
                        'RAW_COUNTS_UNCOLLAPSED', 'SAMPLE_META', 'CELL_SET_AND_POOL_META', 'CELL_LINE_META', 'CONTROL_BARCODE_META',

                        // sushi output files
                        'PRISM_BARCODE_COUNTS', 'UNKNOWN_BARCODE_COUNTS', 'ANNOTATED_COUNTS', 'FILTERED_COUNTS', 'NORMALIZED_COUNTS',
                        'LFC', 'COLLAPSED_LFC',

                        // collate_fastq_reads parameters
                        'SEQUENCING_INDEX_COLS', 'ID_COLS', 'BARCODE_COL', 'LOW_ABUNDANCE_THRESHOLD', 'CHUNK_SIZE', 'REVERSE_INDEX2',

                        // normalize parameters
                        'NORM_METHOD', 'PSEUDOCOUNT', 'READ_DETECTION_LIMIT',

                        // compute_l2fc parameters
                        'SIG_COLS', 'CONTROL_COLS', 'CELL_LINE_COLS', 'COUNT_COL_NAME', 'CTL_TYPES', 'COUNT_THRESHOLD', 'VIABILITY_CAP',
                        'GROWTH_PATTERN_COL', 'BIO_REP_COL',

                        // biomarker parameters
                        'UNIVARIATE_BIOMARKER', 'MULTIVARIATE_BIOMARKER', 'BIOMARKER_FILE', 'DR_COLUMN', 'LFC_BIOMARKER', 'AUC_BIOMARKER',
                        'DR_PATH',

                        //io parameters
                        'MERGE_PATTERNS',

                        //qc parameters
                        'QC_PARAMS', 'FRACTION_EXPECTED_CONTROLS', 'FILTER_FAILED_LINES',

                        //combinations
                        'COMBINATION_COL', 'N_SAMPLES', 'SYNERGY_BIOMARKER', 'SYNERGY_COL', 'SYNERGY_FILE', 'SYNERGY_L2FC_COL'
                    ]

                    def config = [:]

                    // Load existing config if it exists
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    // Add parameters from Jenkins UI only if they don't exist in the config
                    paramList.each { paramName ->
                        if (!config.containsKey(paramName) && params.containsKey(paramName)) {
                            config[paramName] = params[paramName]
                        }
                    }

                    // Add or overwrite specific keys
                    if (!config.containsKey('API_KEY')) {
                        config.API_KEY = sh(script: 'cat /local/jenkins/.clue_api_key', returnStdout: true).trim()
                    }

                    // Explicit settings that are programmatically derived
                    config.REVERSE_INDEX2 = config.SEQ_TYPE == 'DRAGEN'

                    // Write the config back to file after all updates
                    writeFile file: env.CONFIG_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                    echo "Generated config.json: ${config}"
                }
            }
        }

        stage('Generate qc_flag_params') {
            steps {
                script {
                    def paramList = [
                        'nc_variability_threshold', 'error_rate_threshold', 'pc_viability_threshold',
                        'nc_raw_count_threshold', 'contamination_threshold', 'cb_mae_threshold',
                        'cb_spearman_threshold', 'cb_cl_ratio_low_negcon', 'cb_cl_ratio_high_negcon',
                        'cb_cl_ratio_low_poscon', 'cb_cl_ratio_high_poscon', 'well_reads_threshold',
                        'pool_well_delta_threshold', 'pool_well_fraction_threshold', 'fraction_expected_controls'
                    ]

                    def config = [:]

                    // Load existing params if they exist
                    if (fileExists(env.QC_PARAMS_FILE_PATH)) {
                        def configText = readFile(file: env.QC_PARAMS_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    // Add parameters from Jenkins UI only if they don't exist in the config
                    paramList.each { paramName ->
                        if (!config.containsKey(paramName) && params.containsKey(paramName)) {
                            config[paramName] = params[paramName]
                        }
                    }
                    // Write the config back to file after all updates
                    writeFile file: env.QC_PARAMS_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                }
            }
        }

        stage('Add Commit ID to Config') {
            steps {
                script {
                    def commitID = sh(script: 'git rev-parse --short HEAD', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg commit "$commitID" '. + {COMMIT: \$commit}' $configFilePath > ${configFilePath}.tmp && mv ${configFilePath}.tmp $configFilePath
                    """
                    echo "Added commit ID to config.json: $commitID"
                }
            }
        }

        stage('Add Timestamp to Config') {
            steps {
                script {
                    def buildtime = sh(script: 'date -u +"%Y-%m-%dT%H:%M:%SZ"', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg buildtime "$buildtime" '. + {TIMESTAMP: \$buildtime}' $configFilePath > ${configFilePath}.tmp && mv ${configFilePath}.tmp $configFilePath
                    """
                    echo "Added build timestamp to config.json: $buildtime"
                }
            }
        }

        stage('Show available podman containers') {
            steps {
                script {
                    sh '/usr/bin/podman images'
                }
            }
        }

        stage('Run Scripts in Container') {
            steps {
                script {
                    if (params.TRIGGER_BUILD) {
                        def scriptsToRun = []
                        if (params.CREATE_SAMPLE_META) {
                            scriptsToRun.add('create_sample_meta/create_sample_meta.sh')
                        }
                        if (params.CREATE_CELLDB_METADATA) {
                            scriptsToRun.add('create_cell_meta/create_cell_meta.sh')
                        }
                        if (params.COLLATE_FASTQ_READS) {
                            scriptsToRun.add('collate_counts/collate_counts.sh')
                        }
                        if (params.FILTER_COUNTS) {
                            scriptsToRun.add('filter_counts/filter_counts.sh')
                        }
                        if (params.CBNORMALIZE) {
                            scriptsToRun.add('normalize/normalize.sh')
                        }
                        if (params.GENERATE_QC_TABLES) {
                            scriptsToRun.add('qc_tables/qc_tables.sh')
                        }
                        if (params.COMPUTE_LFC) {
                            scriptsToRun.add('compute_l2fc/compute_l2fc.sh')
                        }
                        if (params.BIAS_CORRECTION) {
                            scriptsToRun.add('bias_correction/bias_correction.sh')
                        }
                        if (params.COLLAPSE) {
                            scriptsToRun.add('collapse_replicates/collapse_replicates.sh')
                        }
                        if (params.GENERATE_QC_TABLES_2) {
                            scriptsToRun.add('qc_tables_2/qc_tables_2.sh')
                        }
                        if (params.DRC) {
                            scriptsToRun.add('drc/dose_response.sh')
                        }
                        if (params.SYNERGY) {
                            scriptsToRun.add('synergy/synergy.sh')
                        }
                        if (params.UNIVARIATE_BIOMARKER || params.MULTIVARIATE_BIOMARKER) {
                            scriptsToRun.add('biomarker/biomarker.sh')
                        }
                        if (params.QC_IMAGES) {
                            scriptsToRun.add('filter_counts_qc/filter_counts_qc.sh')
                        }
                        if (params.JOIN_METADATA) {
                            scriptsToRun.add('join_metadata/join_metadata.sh')
                        }
                        if (params.CONVERT_SUSHI) {
                            scriptsToRun.add('sushi_2_s3/sushi_2_s3.sh')
                        }
                        if (params.UPLOAD_TO_BQ) {
                            scriptsToRun.add('load_bq/load_bq.sh')
                        }

                        scriptsToRun.each { scriptName ->
                            echo "Running script: ${scriptName}"

                            sh """
                                chmod +x $WORKSPACE/scripts/launch_job.sh
                                $WORKSPACE/scripts/launch_job.sh $scriptName
                            """
                        }
                    }
                }
            }
        }
    }

    post {
        always {
            script {
                def buildDir = params.BUILD_DIR
                def buildName = params.BUILD_NAME ?: "default_build_name"
                def logFilePath = "${buildDir}/logs/console_output.log"

                // Capture the entire console log
                def log = currentBuild.rawBuild.getLog(999999)

                // Write the log to the specified file
                writeFile file: logFilePath, text: log.join("\n")

                echo "Console output has been saved to ${logFilePath}"
            }
        }
        success {
            script {
                echo 'Build completed successfully.'
            }
        }
    }
}
