#!/usr/bin/env node
import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ErrorCode,
  ListResourcesRequestSchema,
  ListResourceTemplatesRequestSchema,
  ListToolsRequestSchema,
  McpError,
  ReadResourceRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import axios, { AxiosInstance } from 'axios';

// PubChem API interfaces
interface CompoundSearchResult {
  IdentifierList: {
    CID: number[];
  };
}

interface PropertyData {
  PropertyTable: {
    Properties: Array<{
      CID: number;
      MolecularFormula?: string;
      MolecularWeight?: number;
      CanonicalSMILES?: string;
      IsomericSMILES?: string;
      InChI?: string;
      InChIKey?: string;
      IUPACName?: string;
      XLogP?: number;
      TPSA?: number;
      Complexity?: number;
      Charge?: number;
      HBondDonorCount?: number;
      HBondAcceptorCount?: number;
      RotatableBondCount?: number;
      HeavyAtomCount?: number;
      AtomStereoCount?: number;
      DefinedAtomStereoCount?: number;
      BondStereoCount?: number;
      DefinedBondStereoCount?: number;
      Volume3D?: number;
      ConformerCount3D?: number;
    }>;
  };
}

// Type guards and validation functions
const isValidCompoundSearchArgs = (
  args: any
): args is { query: string; search_type?: string; max_records?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.query === 'string' &&
    args.query.length > 0 &&
    (args.search_type === undefined || ['name', 'smiles', 'inchi', 'sdf', 'cid', 'formula'].includes(args.search_type)) &&
    (args.max_records === undefined || (typeof args.max_records === 'number' && args.max_records > 0 && args.max_records <= 10000))
  );
};

const isValidCidArgs = (
  args: any
): args is { cid: number | string; format?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.format === undefined || ['json', 'sdf', 'xml', 'asnt', 'asnb'].includes(args.format))
  );
};

const isValidSmilesArgs = (
  args: any
): args is { smiles: string; threshold?: number; max_records?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.smiles === 'string' &&
    args.smiles.length > 0 &&
    (args.threshold === undefined || (typeof args.threshold === 'number' && args.threshold >= 0 && args.threshold <= 100)) &&
    (args.max_records === undefined || (typeof args.max_records === 'number' && args.max_records > 0 && args.max_records <= 10000))
  );
};

const isValidBatchArgs = (
  args: any
): args is { cids: number[]; operation?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    Array.isArray(args.cids) &&
    args.cids.length > 0 &&
    args.cids.length <= 200 &&
    args.cids.every((cid: any) => typeof cid === 'number' && cid > 0) &&
    (args.operation === undefined || ['property', 'synonyms', 'classification', 'description'].includes(args.operation))
  );
};

const isValidConformerArgs = (
  args: any
): args is { cid: number | string; conformer_type?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.conformer_type === undefined || ['3d', '2d'].includes(args.conformer_type))
  );
};

const isValidPropertiesArgs = (
  args: any
): args is { cid: number | string; properties?: string[] } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.properties === undefined || (Array.isArray(args.properties) && args.properties.every((p: any) => typeof p === 'string')))
  );
};

class PubChemServer {
  private server: Server;
  private apiClient: AxiosInstance;

  constructor() {
    this.server = new Server(
      {
        name: 'pubchem-server',
        version: '1.0.0',
      },
      {
        capabilities: {
          resources: {},
          tools: {},
        },
      }
    );

    // Initialize PubChem API client
    this.apiClient = axios.create({
      baseURL: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug',
      timeout: 30000,
      headers: {
        'User-Agent': 'PubChem-MCP-Server/1.0.0',
        'Accept': 'application/json',
      },
    });

    this.setupResourceHandlers();
    this.setupToolHandlers();

    // Error handling
    this.server.onerror = (error: any) => console.error('[MCP Error]', error);
    process.on('SIGINT', async () => {
      await this.server.close();
      process.exit(0);
    });
  }

  private setupResourceHandlers() {
    // List available resource templates
    this.server.setRequestHandler(
      ListResourceTemplatesRequestSchema,
      async () => ({
        resourceTemplates: [
          {
            uriTemplate: 'pubchem://compound/{cid}',
            name: 'PubChem compound entry',
            mimeType: 'application/json',
            description: 'Complete compound information for a PubChem CID',
          },
          {
            uriTemplate: 'pubchem://structure/{cid}',
            name: 'Chemical structure data',
            mimeType: 'application/json',
            description: '2D/3D structure information for a compound',
          },
          {
            uriTemplate: 'pubchem://properties/{cid}',
            name: 'Chemical properties',
            mimeType: 'application/json',
            description: 'Molecular properties and descriptors for a compound',
          },
          {
            uriTemplate: 'pubchem://bioassay/{aid}',
            name: 'PubChem bioassay data',
            mimeType: 'application/json',
            description: 'Bioassay information and results for an AID',
          },
          {
            uriTemplate: 'pubchem://similarity/{smiles}',
            name: 'Similarity search results',
            mimeType: 'application/json',
            description: 'Chemical similarity search results for a SMILES string',
          },
          {
            uriTemplate: 'pubchem://safety/{cid}',
            name: 'Safety and toxicity data',
            mimeType: 'application/json',
            description: 'Safety classifications and toxicity information',
          },
        ],
      })
    );

    // Handle resource requests
    this.server.setRequestHandler(
      ReadResourceRequestSchema,
      async (request: any) => {
        const uri = request.params.uri;

        // Handle compound info requests
        const compoundMatch = uri.match(/^pubchem:\/\/compound\/([0-9]+)$/);
        if (compoundMatch) {
          const cid = compoundMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch compound ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle structure requests
        const structureMatch = uri.match(/^pubchem:\/\/structure\/([0-9]+)$/);
        if (structureMatch) {
          const cid = structureMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch structure for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle properties requests
        const propertiesMatch = uri.match(/^pubchem:\/\/properties\/([0-9]+)$/);
        if (propertiesMatch) {
          const cid = propertiesMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,Complexity/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch properties for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle bioassay requests
        const bioassayMatch = uri.match(/^pubchem:\/\/bioassay\/([0-9]+)$/);
        if (bioassayMatch) {
          const aid = bioassayMatch[1];
          try {
            const response = await this.apiClient.get(`/assay/aid/${aid}/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch bioassay ${aid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle similarity search requests
        const similarityMatch = uri.match(/^pubchem:\/\/similarity\/(.+)$/);
        if (similarityMatch) {
          const smiles = decodeURIComponent(similarityMatch[1]);
          try {
            const response = await this.apiClient.post('/compound/similarity/smiles/JSON', {
              smiles: smiles,
              Threshold: 90,
              MaxRecords: 100,
            });
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to perform similarity search: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle safety data requests
        const safetyMatch = uri.match(/^pubchem:\/\/safety\/([0-9]+)$/);
        if (safetyMatch) {
          const cid = safetyMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/classification/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch safety data for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        throw new McpError(
          ErrorCode.InvalidRequest,
          `Invalid URI format: ${uri}`
        );
      }
    );
  }

  private setupToolHandlers() {
    this.server.setRequestHandler(ListToolsRequestSchema, async () => ({
      tools: [
        // Chemical Search & Retrieval (6 tools)
        {
          name: 'search_compounds',
          description: 'Search PubChem database for compounds by name, CAS number, formula, or identifier',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'Search query (compound name, CAS, formula, or identifier)' },
              search_type: { type: 'string', enum: ['name', 'smiles', 'inchi', 'sdf', 'cid', 'formula'], description: 'Type of search to perform (default: name)' },
              max_records: { type: 'number', description: 'Maximum number of results (1-10000, default: 100)', minimum: 1, maximum: 10000 },
            },
            required: ['query'],
          },
        },
        {
          name: 'get_compound_info',
          description: 'Get detailed information for a specific compound by PubChem CID',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              format: { type: 'string', enum: ['json', 'sdf', 'xml', 'asnt', 'asnb'], description: 'Output format (default: json)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'search_by_smiles',
          description: 'Search for compounds by SMILES string (exact match)',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the query molecule' },
            },
            required: ['smiles'],
          },
        },
        {
          name: 'search_by_inchi',
          description: 'Search for compounds by InChI or InChI key',
          inputSchema: {
            type: 'object',
            properties: {
              inchi: { type: 'string', description: 'InChI string or InChI key' },
            },
            required: ['inchi'],
          },
        },
        {
          name: 'search_by_cas_number',
          description: 'Search for compounds by CAS Registry Number',
          inputSchema: {
            type: 'object',
            properties: {
              cas_number: { type: 'string', description: 'CAS Registry Number (e.g., 50-78-2)' },
            },
            required: ['cas_number'],
          },
        },
        {
          name: 'get_compound_synonyms',
          description: 'Get all names and synonyms for a compound',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },

        // Structure Analysis & Similarity (5 tools)
        {
          name: 'search_similar_compounds',
          description: 'Find chemically similar compounds using Tanimoto similarity',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the query molecule' },
              threshold: { type: 'number', description: 'Similarity threshold (0-100, default: 90)', minimum: 0, maximum: 100 },
              max_records: { type: 'number', description: 'Maximum number of results (1-10000, default: 100)', minimum: 1, maximum: 10000 },
            },
            required: ['smiles'],
          },
        },
        {
          name: 'substructure_search',
          description: 'Find compounds containing a specific substructure',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the substructure query' },
              max_records: { type: 'number', description: 'Maximum number of results (1-10000, default: 100)', minimum: 1, maximum: 10000 },
            },
            required: ['smiles'],
          },
        },
        {
          name: 'superstructure_search',
          description: 'Find larger compounds that contain the query structure',
          inputSchema: {
            type: 'object',
            properties: {
              smiles: { type: 'string', description: 'SMILES string of the query structure' },
              max_records: { type: 'number', description: 'Maximum number of results (1-10000, default: 100)', minimum: 1, maximum: 10000 },
            },
            required: ['smiles'],
          },
        },
        {
          name: 'get_3d_conformers',
          description: 'Get 3D conformer data and structural information',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              conformer_type: { type: 'string', enum: ['3d', '2d'], description: 'Type of conformer data (default: 3d)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'analyze_stereochemistry',
          description: 'Analyze stereochemistry, chirality, and isomer information',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },

        // Chemical Properties & Descriptors (6 tools)
        {
          name: 'get_compound_properties',
          description: 'Get molecular properties (MW, logP, TPSA, etc.)',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              properties: { type: 'array', items: { type: 'string' }, description: 'Specific properties to retrieve (optional)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'calculate_descriptors',
          description: 'Calculate comprehensive molecular descriptors and fingerprints',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              descriptor_type: { type: 'string', enum: ['all', 'basic', 'topological', '3d'], description: 'Type of descriptors (default: all)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'predict_admet_properties',
          description: 'Predict ADMET properties (Absorption, Distribution, Metabolism, Excretion, Toxicity)',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              smiles: { type: 'string', description: 'SMILES string (alternative to CID)' },
            },
            required: [],
          },
        },
        {
          name: 'assess_drug_likeness',
          description: 'Assess drug-likeness using Lipinski Rule of Five, Veber rules, and PAINS filters',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              smiles: { type: 'string', description: 'SMILES string (alternative to CID)' },
            },
            required: [],
          },
        },
        {
          name: 'analyze_molecular_complexity',
          description: 'Analyze molecular complexity and synthetic accessibility',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'get_pharmacophore_features',
          description: 'Get pharmacophore features and binding site information',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },

        // Bioassay & Activity Data (5 tools)
        {
          name: 'search_bioassays',
          description: 'Search for biological assays by target, description, or source',
          inputSchema: {
            type: 'object',
            properties: {
              query: { type: 'string', description: 'General search query' },
              target: { type: 'string', description: 'Target protein or gene name' },
              source: { type: 'string', description: 'Data source (e.g., ChEMBL, NCGC)' },
              max_records: { type: 'number', description: 'Maximum number of results (1-1000, default: 100)', minimum: 1, maximum: 1000 },
            },
            required: [],
          },
        },
        {
          name: 'get_assay_info',
          description: 'Get detailed information for a specific bioassay by AID',
          inputSchema: {
            type: 'object',
            properties: {
              aid: { type: 'number', description: 'PubChem Assay ID (AID)' },
            },
            required: ['aid'],
          },
        },
        {
          name: 'get_compound_bioactivities',
          description: 'Get all bioassay results and activities for a compound',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              activity_outcome: { type: 'string', enum: ['active', 'inactive', 'inconclusive', 'all'], description: 'Filter by activity outcome (default: all)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'search_by_target',
          description: 'Find compounds tested against a specific biological target',
          inputSchema: {
            type: 'object',
            properties: {
              target: { type: 'string', description: 'Target name (gene, protein, or pathway)' },
              activity_type: { type: 'string', description: 'Type of activity (e.g., IC50, EC50, Ki)' },
              max_records: { type: 'number', description: 'Maximum number of results (1-1000, default: 100)', minimum: 1, maximum: 1000 },
            },
            required: ['target'],
          },
        },
        {
          name: 'compare_activity_profiles',
          description: 'Compare bioactivity profiles across multiple compounds',
          inputSchema: {
            type: 'object',
            properties: {
              cids: { type: 'array', items: { type: 'number' }, description: 'Array of PubChem CIDs (2-50)', minItems: 2, maxItems: 50 },
              activity_type: { type: 'string', description: 'Specific activity type for comparison (optional)' },
            },
            required: ['cids'],
          },
        },

        // Safety & Toxicity (4 tools)
        {
          name: 'get_safety_data',
          description: 'Get GHS hazard classifications and safety information',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'get_toxicity_info',
          description: 'Get toxicity data including LD50, carcinogenicity, and mutagenicity',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'assess_environmental_fate',
          description: 'Assess environmental fate including biodegradation and bioaccumulation',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'get_regulatory_info',
          description: 'Get regulatory information from FDA, EPA, and international agencies',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },

        // Cross-References & Integration (4 tools)
        {
          name: 'get_external_references',
          description: 'Get links to external databases (ChEMBL, DrugBank, KEGG, etc.)',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'search_patents',
          description: 'Search for chemical patents and intellectual property information',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
              query: { type: 'string', description: 'Patent search query (alternative to CID)' },
            },
            required: [],
          },
        },
        {
          name: 'get_literature_references',
          description: 'Get PubMed citations and scientific literature references',
          inputSchema: {
            type: 'object',
            properties: {
              cid: { type: ['number', 'string'], description: 'PubChem Compound ID (CID)' },
            },
            required: ['cid'],
          },
        },
        {
          name: 'batch_compound_lookup',
          description: 'Process multiple compound IDs efficiently',
          inputSchema: {
            type: 'object',
            properties: {
              cids: { type: 'array', items: { type: 'number' }, description: 'Array of PubChem CIDs (1-200)', minItems: 1, maxItems: 200 },
              operation: { type: 'string', enum: ['property', 'synonyms', 'classification', 'description'], description: 'Operation to perform (default: property)' },
            },
            required: ['cids'],
          },
        },
      ],
    }));

    this.server.setRequestHandler(CallToolRequestSchema, async (request: any) => {
      const { name, arguments: args } = request.params;

      try {
        switch (name) {
          // Chemical Search & Retrieval
          case 'search_compounds':
            return await this.handleSearchCompounds(args);
          case 'get_compound_info':
            return await this.handleGetCompoundInfo(args);
          case 'search_by_smiles':
            return await this.handleSearchBySmiles(args);
          case 'search_by_inchi':
            return await this.handleSearchByInchi(args);
          case 'search_by_cas_number':
            return await this.handleSearchByCasNumber(args);
          case 'get_compound_synonyms':
            return await this.handleGetCompoundSynonyms(args);

          // Structure Analysis & Similarity
          case 'search_similar_compounds':
            return await this.handleSearchSimilarCompounds(args);
          case 'substructure_search':
            return await this.handleSubstructureSearch(args);
          case 'superstructure_search':
            return await this.handleSuperstructureSearch(args);
          case 'get_3d_conformers':
            return await this.handleGet3dConformers(args);
          case 'analyze_stereochemistry':
            return await this.handleAnalyzeStereochemistry(args);

          // Chemical Properties & Descriptors
          case 'get_compound_properties':
            return await this.handleGetCompoundProperties(args);
          case 'calculate_descriptors':
            return await this.handleCalculateDescriptors(args);
          case 'predict_admet_properties':
            return await this.handlePredictAdmetProperties(args);
          case 'assess_drug_likeness':
            return await this.handleAssessDrugLikeness(args);
          case 'analyze_molecular_complexity':
            return await this.handleAnalyzeMolecularComplexity(args);
          case 'get_pharmacophore_features':
            return await this.handleGetPharmacophoreFeatures(args);

          // Bioassay & Activity Data
          case 'search_bioassays':
            return await this.handleSearchBioassays(args);
          case 'get_assay_info':
            return await this.handleGetAssayInfo(args);
          case 'get_compound_bioactivities':
            return await this.handleGetCompoundBioactivities(args);
          case 'search_by_target':
            return await this.handleSearchByTarget(args);
          case 'compare_activity_profiles':
            return await this.handleCompareActivityProfiles(args);

          // Safety & Toxicity
          case 'get_safety_data':
            return await this.handleGetSafetyData(args);
          case 'get_toxicity_info':
            return await this.handleGetToxicityInfo(args);
          case 'assess_environmental_fate':
            return await this.handleAssessEnvironmentalFate(args);
          case 'get_regulatory_info':
            return await this.handleGetRegulatoryInfo(args);

          // Cross-References & Integration
          case 'get_external_references':
            return await this.handleGetExternalReferences(args);
          case 'search_patents':
            return await this.handleSearchPatents(args);
          case 'get_literature_references':
            return await this.handleGetLiteratureReferences(args);
          case 'batch_compound_lookup':
            return await this.handleBatchCompoundLookup(args);

          default:
            throw new McpError(
              ErrorCode.MethodNotFound,
              `Unknown tool: ${name}`
            );
        }
      } catch (error) {
        return {
          content: [
            {
              type: 'text',
              text: `Error executing tool ${name}: ${error instanceof Error ? error.message : 'Unknown error'}`,
            },
          ],
          isError: true,
        };
      }
    });
  }

  // Chemical Search & Retrieval handlers
  private async handleSearchCompounds(args: any) {
    if (!isValidCompoundSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid compound search arguments');
    }

    try {
      const searchType = args.search_type || 'name';
      const maxRecords = args.max_records || 100;

      const response = await this.apiClient.get(`/compound/${searchType}/${encodeURIComponent(args.query)}/cids/JSON`, {
        params: {
          MaxRecords: maxRecords,
        },
      });

      if (response.data?.IdentifierList?.CID?.length > 0) {
        const cids = response.data.IdentifierList.CID.slice(0, 10);
        const detailsResponse = await this.apiClient.get(`/compound/cid/${cids.join(',')}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                query: args.query,
                search_type: searchType,
                total_found: response.data.IdentifierList.CID.length,
                details: detailsResponse.data,
              }, null, 2),
            },
          ],
        };
      }

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({ message: 'No compounds found', query: args.query }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundInfo(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const format = args.format || 'json';
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/${format === 'json' ? 'JSON' : format}`);

      return {
        content: [
          {
            type: 'text',
            text: format === 'json' ? JSON.stringify(response.data, null, 2) : String(response.data),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound info: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchBySmiles(args: any) {
    if (!isValidSmilesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid SMILES arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/smiles/${encodeURIComponent(args.smiles)}/cids/JSON`);

      if (response.data?.IdentifierList?.CID?.length > 0) {
        const cid = response.data.IdentifierList.CID[0];
        const detailsResponse = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                query_smiles: args.smiles,
                found_cid: cid,
                details: detailsResponse.data,
              }, null, 2),
            },
          ],
        };
      }

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({ message: 'No exact match found', query_smiles: args.smiles }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search by SMILES: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // Simplified implementation handlers (placeholder implementations)
  private async handleSearchByInchi(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'InChI search not yet implemented', args }, null, 2) }] };
  }

  private async handleSearchByCasNumber(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'CAS search not yet implemented', args }, null, 2) }] };
  }

  private async handleGetCompoundSynonyms(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/synonyms/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound synonyms: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchSimilarCompounds(args: any) {
    if (!isValidSmilesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid similarity search arguments');
    }

    try {
      const threshold = args.threshold || 90;
      const maxRecords = args.max_records || 100;

      const response = await this.apiClient.post('/compound/similarity/smiles/JSON', {
        smiles: args.smiles,
        Threshold: threshold,
        MaxRecords: maxRecords,
      });

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search similar compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSubstructureSearch(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Substructure search not yet implemented', args }, null, 2) }] };
  }

  private async handleSuperstructureSearch(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Superstructure search not yet implemented', args }, null, 2) }] };
  }

  private async handleGet3dConformers(args: any) {
    if (!isValidConformerArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid 3D conformer arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/Volume3D,ConformerCount3D/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              cid: args.cid,
              conformer_type: args.conformer_type || '3d',
              properties: response.data,
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get 3D conformers: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleAnalyzeStereochemistry(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid stereochemistry arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/AtomStereoCount,DefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,IsomericSMILES/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              cid: args.cid,
              stereochemistry: response.data,
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to analyze stereochemistry: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundProperties(args: any) {
    if (!isValidPropertiesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid compound properties arguments');
    }

    try {
      const properties = args.properties || [
        'MolecularWeight', 'XLogP', 'TPSA', 'HBondDonorCount', 'HBondAcceptorCount',
        'RotatableBondCount', 'Complexity', 'HeavyAtomCount', 'Charge'
      ];

      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/${properties.join(',')}/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound properties: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // Placeholder implementations for remaining methods
  private async handleCalculateDescriptors(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Descriptor calculation not yet implemented', args }, null, 2) }] };
  }

  private async handlePredictAdmetProperties(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'ADMET prediction not yet implemented', args }, null, 2) }] };
  }

  private async handleAssessDrugLikeness(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Drug-likeness assessment not yet implemented', args }, null, 2) }] };
  }

  private async handleAnalyzeMolecularComplexity(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Molecular complexity analysis not yet implemented', args }, null, 2) }] };
  }

  private async handleGetPharmacophoreFeatures(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Pharmacophore features not yet implemented', args }, null, 2) }] };
  }

  private async handleSearchBioassays(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Bioassay search not yet implemented', args }, null, 2) }] };
  }

  private async handleGetAssayInfo(args: any) {
    try {
      const response = await this.apiClient.get(`/assay/aid/${args.aid}/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get assay info: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundBioactivities(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Bioactivity search not yet implemented', args }, null, 2) }] };
  }

  private async handleSearchByTarget(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Target search not yet implemented', args }, null, 2) }] };
  }

  private async handleCompareActivityProfiles(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Activity profile comparison not yet implemented', args }, null, 2) }] };
  }

  private async handleGetSafetyData(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/classification/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get safety data: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetToxicityInfo(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const cid = args.cid;
      const toxicityData: any = {
        cid: cid,
        toxicity_summary: {
          acute_toxicity: {},
          chronic_toxicity: {},
          carcinogenicity: {},
          mutagenicity: {},
          reproductive_toxicity: {},
          environmental_toxicity: {}
        },
        data_sources: [],
        last_updated: new Date().toISOString()
      };

      // Try to get GHS classification data which includes toxicity information
      try {
        const ghsResponse = await this.apiClient.get(`/compound/cid/${cid}/classification/JSON`);
        if (ghsResponse.data?.Informations) {
          const ghsInfo = ghsResponse.data.Informations[0];
          if (ghsInfo?.GHSClassification) {
            toxicityData.ghs_classification = ghsInfo.GHSClassification;
            toxicityData.data_sources.push('GHS Classification');
          }
        }
      } catch (error) {
        // GHS data not available, continue with other sources
      }

      // Try to get compound description which may contain toxicity information
      try {
        const descResponse = await this.apiClient.get(`/compound/cid/${cid}/description/JSON`);
        if (descResponse.data?.InformationList?.Information) {
          const descriptions = descResponse.data.InformationList.Information;
          const toxicityDescriptions = descriptions.filter((desc: any) => 
            desc.Title && (
              desc.Title.toLowerCase().includes('toxicity') ||
              desc.Title.toLowerCase().includes('toxic') ||
              desc.Title.toLowerCase().includes('ld50') ||
              desc.Title.toLowerCase().includes('carcinogen') ||
              desc.Title.toLowerCase().includes('mutagen')
            )
          );
          
          if (toxicityDescriptions.length > 0) {
            toxicityData.toxicity_descriptions = toxicityDescriptions.map((desc: any) => ({
              title: desc.Title,
              description: desc.Description,
              source: desc.SourceName,
              url: desc.ReferenceNumber ? `https://pubchem.ncbi.nlm.nih.gov/source/${desc.ReferenceNumber}` : null
            }));
            toxicityData.data_sources.push('PubChem Descriptions');
          }
        }
      } catch (error) {
        // Description data not available
      }

      // Get basic compound properties that relate to toxicity
      try {
        const propResponse = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON`);
        if (propResponse.data?.PropertyTable?.Properties?.[0]) {
          const props = propResponse.data.PropertyTable.Properties[0];
          toxicityData.molecular_properties = {
            molecular_weight: props.MolecularWeight,
            logp: props.XLogP,
            polar_surface_area: props.TPSA,
            hbd_count: props.HBondDonorCount,
            hba_count: props.HBondAcceptorCount
          };
          
          // Basic toxicity predictions based on molecular properties
          toxicityData.predicted_toxicity = {
            oral_toxicity_class: this.predictOralToxicity(props.MolecularWeight, props.XLogP),
            bioaccumulation_potential: props.XLogP > 3 ? 'High' : props.XLogP > 1 ? 'Moderate' : 'Low',
            membrane_permeability: props.TPSA < 140 ? 'High' : 'Low'
          };
        }
      } catch (error) {
        // Properties not available
      }

      // Add general toxicity assessment notes
      toxicityData.assessment_notes = [
        'This is a computational assessment based on available PubChem data',
        'For definitive toxicity information, consult official regulatory databases',
        'Experimental validation is required for safety assessments',
        'Data availability varies significantly between compounds'
      ];

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(toxicityData, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error fetching toxicity info: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleAssessEnvironmentalFate(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const cid = args.cid;
      const environmentalData: any = {
        cid: cid,
        environmental_fate: {
          biodegradation: {},
          bioaccumulation: {},
          atmospheric_fate: {},
          aquatic_fate: {},
          soil_fate: {}
        },
        regulatory_status: {},
        data_sources: [],
        last_updated: new Date().toISOString()
      };

      // Get molecular properties for environmental fate prediction
      try {
        const propResponse = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,XLogP,TPSA,CanonicalSMILES,MolecularFormula/JSON`);
        if (propResponse.data?.PropertyTable?.Properties?.[0]) {
          const props = propResponse.data.PropertyTable.Properties[0];
          environmentalData.molecular_properties = {
            molecular_weight: props.MolecularWeight,
            logp: props.XLogP,
            polar_surface_area: props.TPSA,
            smiles: props.CanonicalSMILES,
            formula: props.MolecularFormula
          };

          // Predict environmental fate based on molecular properties
          environmentalData.environmental_fate = {
            biodegradation: {
              predicted_rate: this.predictBiodegradation(props.MolecularWeight, props.XLogP),
              factors: [
                'Molecular weight affects microbial uptake',
                'LogP influences bioavailability',
                'Structural complexity affects degradation pathways'
              ]
            },
            bioaccumulation: {
              potential: props.XLogP > 3 ? 'High' : props.XLogP > 1 ? 'Moderate' : 'Low',
              bcf_estimate: props.XLogP ? Math.pow(10, 0.85 * props.XLogP - 0.70) : 'Not calculated',
              factors: ['LogP is primary indicator of bioaccumulation potential']
            },
            atmospheric_fate: {
              volatility: props.MolecularWeight < 200 ? 'High' : props.MolecularWeight < 500 ? 'Moderate' : 'Low',
              photodegradation: 'Depends on molecular structure and UV absorption',
              factors: ['Molecular weight affects vapor pressure']
            },
            aquatic_fate: {
              solubility_prediction: props.XLogP < 0 ? 'High' : props.XLogP < 3 ? 'Moderate' : 'Low',
              sorption_potential: props.XLogP > 2.5 ? 'High' : 'Low',
              factors: ['LogP correlates with water solubility and sorption']
            },
            soil_fate: {
              mobility: props.XLogP < 2 ? 'High' : props.XLogP < 4 ? 'Moderate' : 'Low',
              persistence: 'Depends on soil conditions and microbial activity',
              factors: ['Sorption to organic matter increases with LogP']
            }
          };

          environmentalData.data_sources.push('Molecular Property Predictions');
        }
      } catch (error) {
        // Properties not available
      }

      // Try to get environmental classification data
      try {
        const classResponse = await this.apiClient.get(`/compound/cid/${cid}/classification/JSON`);
        if (classResponse.data?.Informations) {
          const classInfo = classResponse.data.Informations[0];
          if (classInfo?.GHSClassification) {
            // Look for environmental hazard classifications
            const envHazards = classInfo.GHSClassification.filter((item: any) => 
              item.Category && (
                item.Category.toLowerCase().includes('aquatic') ||
                item.Category.toLowerCase().includes('environmental') ||
                item.Category.toLowerCase().includes('ozone')
              )
            );
            
            if (envHazards.length > 0) {
              environmentalData.environmental_hazards = envHazards;
              environmentalData.data_sources.push('GHS Environmental Classification');
            }
          }
        }
      } catch (error) {
        // Classification data not available
      }

      // Add assessment methodology and limitations
      environmentalData.assessment_methodology = {
        approach: 'QSAR-based predictions using molecular descriptors',
        limitations: [
          'Predictions are estimates based on molecular properties',
          'Actual environmental fate depends on specific conditions',
          'Experimental data should be used when available',
          'Metabolites and transformation products not considered'
        ],
        recommended_tests: [
          'OECD biodegradation studies',
          'Bioaccumulation factor determination',
          'Environmental monitoring data',
          'Ecotoxicity testing'
        ]
      };

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(environmentalData, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error assessing environmental fate: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  private async handleGetRegulatoryInfo(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const cid = args.cid;
      const regulatoryData: any = {
        cid: cid,
        regulatory_status: {
          fda: {},
          epa: {},
          echa: {},
          international: {}
        },
        classifications: {},
        data_sources: [],
        last_updated: new Date().toISOString()
      };

      // Get basic compound information
      try {
        const compoundResponse = await this.apiClient.get(`/compound/cid/${cid}/JSON`);
        if (compoundResponse.data?.PC_Compounds?.[0]) {
          const compound = compoundResponse.data.PC_Compounds[0];
          regulatoryData.compound_info = {
            iupac_name: compound.props?.find((p: any) => p.urn?.label === 'IUPAC Name')?.value?.sval,
            molecular_formula: compound.props?.find((p: any) => p.urn?.label === 'Molecular Formula')?.value?.sval,
            cas_number: compound.props?.find((p: any) => p.urn?.label === 'CAS')?.value?.sval
          };
        }
      } catch (error) {
        // Compound info not available
      }

      // Get GHS classification which includes regulatory classifications
      try {
        const ghsResponse = await this.apiClient.get(`/compound/cid/${cid}/classification/JSON`);
        if (ghsResponse.data?.Informations) {
          const ghsInfo = ghsResponse.data.Informations[0];
          if (ghsInfo?.GHSClassification) {
            regulatoryData.classifications.ghs = {
              hazard_classes: ghsInfo.GHSClassification.map((item: any) => ({
                category: item.Category,
                classification: item.Classification,
                hazard_statement: item.HazardStatement,
                signal_word: item.SignalWord
              })),
              source: ghsInfo.SourceName,
              reference: ghsInfo.ReferenceNumber
            };
            regulatoryData.data_sources.push('GHS Classification System');
          }
        }
      } catch (error) {
        // GHS data not available
      }

      // Try to get compound descriptions that may contain regulatory information
      try {
        const descResponse = await this.apiClient.get(`/compound/cid/${cid}/description/JSON`);
        if (descResponse.data?.InformationList?.Information) {
          const descriptions = descResponse.data.InformationList.Information;
          const regulatoryDescriptions = descriptions.filter((desc: any) => 
            desc.Title && (
              desc.Title.toLowerCase().includes('regulation') ||
              desc.Title.toLowerCase().includes('regulatory') ||
              desc.Title.toLowerCase().includes('fda') ||
              desc.Title.toLowerCase().includes('epa') ||
              desc.Title.toLowerCase().includes('approval') ||
              desc.Title.toLowerCase().includes('controlled') ||
              desc.Title.toLowerCase().includes('schedule')
            )
          );
          
          if (regulatoryDescriptions.length > 0) {
            regulatoryData.regulatory_descriptions = regulatoryDescriptions.map((desc: any) => ({
              title: desc.Title,
              description: desc.Description,
              source: desc.SourceName,
              url: desc.ReferenceNumber ? `https://pubchem.ncbi.nlm.nih.gov/source/${desc.ReferenceNumber}` : null
            }));
            regulatoryData.data_sources.push('PubChem Regulatory Descriptions');
          }
        }
      } catch (error) {
        // Description data not available
      }

      // Add regulatory framework information
      regulatoryData.regulatory_frameworks = {
        united_states: {
          fda: {
            description: 'Food and Drug Administration',
            scope: 'Food additives, drugs, cosmetics, medical devices',
            database_url: 'https://www.fda.gov/'
          },
          epa: {
            description: 'Environmental Protection Agency',
            scope: 'Pesticides, industrial chemicals, environmental contaminants',
            database_url: 'https://www.epa.gov/'
          },
          osha: {
            description: 'Occupational Safety and Health Administration',
            scope: 'Workplace chemical safety',
            database_url: 'https://www.osha.gov/'
          }
        },
        european_union: {
          echa: {
            description: 'European Chemicals Agency',
            scope: 'REACH regulation, CLP classification',
            database_url: 'https://echa.europa.eu/'
          },
          ema: {
            description: 'European Medicines Agency',
            scope: 'Pharmaceutical regulation',
            database_url: 'https://www.ema.europa.eu/'
          }
        },
        international: {
          who: {
            description: 'World Health Organization',
            scope: 'Global health standards',
            database_url: 'https://www.who.int/'
          },
          oecd: {
            description: 'Organisation for Economic Co-operation and Development',
            scope: 'Chemical testing guidelines',
            database_url: 'https://www.oecd.org/'
          }
        }
      };

      // Add data limitations and recommendations
      regulatoryData.data_limitations = [
        'PubChem aggregates data from multiple sources with varying update frequencies',
        'Regulatory status can change frequently - verify with official sources',
        'Some regulatory information may not be publicly available',
        'Regional regulations may differ significantly'
      ];

      regulatoryData.recommendations = [
        'Consult official regulatory databases for current status',
        'Check multiple jurisdictions for comprehensive coverage',
        'Consider consulting regulatory affairs professionals',
        'Monitor regulatory updates for changes in status'
      ];

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(regulatoryData, null, 2),
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: 'text',
            text: `Error fetching regulatory info: ${error instanceof Error ? error.message : 'Unknown error'}`,
          },
        ],
        isError: true,
      };
    }
  }

  // Helper methods for toxicity and environmental predictions
  private predictOralToxicity(molecularWeight?: number, logP?: number): string {
    if (!molecularWeight || !logP) return 'Cannot predict - insufficient data';
    
    // Simple QSAR-based prediction
    if (molecularWeight > 500 && logP > 5) return 'Class IV (Low toxicity)';
    if (molecularWeight < 200 && logP < 0) return 'Class II (Moderate toxicity)';
    if (logP > 3 && molecularWeight < 300) return 'Class III (Moderate-low toxicity)';
    return 'Class III-IV (Moderate-low toxicity)';
  }

  private predictBiodegradation(molecularWeight?: number, logP?: number): string {
    if (!molecularWeight || !logP) return 'Cannot predict - insufficient data';
    
    if (molecularWeight < 200 && logP < 3) return 'Fast (days to weeks)';
    if (molecularWeight < 500 && logP < 4) return 'Moderate (weeks to months)';
    if (molecularWeight > 500 || logP > 5) return 'Slow (months to years)';
    return 'Moderate (weeks to months)';
  }

  private async handleGetExternalReferences(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'External references not yet implemented', args }, null, 2) }] };
  }

  private async handleSearchPatents(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Patent search not yet implemented', args }, null, 2) }] };
  }

  private async handleGetLiteratureReferences(args: any) {
    return { content: [{ type: 'text', text: JSON.stringify({ message: 'Literature references not yet implemented', args }, null, 2) }] };
  }

  private async handleBatchCompoundLookup(args: any) {
    if (!isValidBatchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid batch arguments');
    }

    try {
      const results = [];
      for (const cid of args.cids.slice(0, 10)) {
        try {
          const response = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,CanonicalSMILES,IUPACName/JSON`);
          results.push({ cid, data: response.data, success: true });
        } catch (error) {
          results.push({ cid, error: error instanceof Error ? error.message : 'Unknown error', success: false });
        }
      }

      return { content: [{ type: 'text', text: JSON.stringify({ batch_results: results }, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Batch lookup failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  async run() {
    const transport = new StdioServerTransport();
    await this.server.connect(transport);
    console.error('PubChem MCP server running on stdio');
  }
}

const server = new PubChemServer();
server.run().catch(console.error);
