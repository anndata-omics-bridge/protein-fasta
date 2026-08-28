"""Backend-owned registry tables, indexes, and relational invariants.

Kept apart from the code that runs it so a backend implementation can read the
DDL without importing the module that opens connections, and so the shape of the
registry can be read in one place.
"""

from __future__ import annotations

STATS_ONLY_TABLES = ("databases", "database_kind_stats", "database_pair_stats", "registry_meta")
"""Registry tables the read paths need; ``databases`` first so references resolve."""

REGISTRY_TABLES: dict[str, str] = {
    "registry_meta": """
        CREATE TABLE IF NOT EXISTS registry_meta (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            ) WITHOUT ROWID
        """,
    "databases": """
        CREATE TABLE IF NOT EXISTS databases (
                id INTEGER PRIMARY KEY,
                relative_path TEXT NOT NULL UNIQUE,
                filename TEXT NOT NULL,
                dbname TEXT NOT NULL,
                file_size_bytes INTEGER NOT NULL,
                mtime_ns INTEGER NOT NULL,
                sentinel_header TEXT,
                annotation TEXT,
                filename_is_decoy INTEGER NOT NULL CHECK (filename_is_decoy IN (0, 1)),
                is_decoy INTEGER NOT NULL CHECK (is_decoy IN (0, 1)),
                contaminant_markers_json TEXT NOT NULL,
                indexed_at TEXT NOT NULL,
                detail_level TEXT NOT NULL CHECK (detail_level IN ('full', 'metadata_only')),
                entry_count INTEGER NOT NULL,
                target_count INTEGER NOT NULL,
                decoy_count INTEGER NOT NULL,
                contaminant_count INTEGER NOT NULL,
                entrapment_count INTEGER NOT NULL,
                sentinel_count INTEGER NOT NULL,
                distinct_target_ids INTEGER,
                distinct_target_sequences INTEGER,
                distinct_target_descriptions INTEGER,
                duplicate_target_id_occurrences INTEGER,
                conflicting_target_ids INTEGER,
                repeated_target_sequences INTEGER,
                length_min INTEGER NOT NULL,
                length_q1 REAL NOT NULL,
                length_median REAL NOT NULL,
                length_mean REAL NOT NULL,
                length_q3 REAL NOT NULL,
                length_max INTEGER NOT NULL,
                total_residues INTEGER NOT NULL,
                aa_sample_size INTEGER NOT NULL CHECK (
                    aa_sample_size >= 0 AND aa_sample_size <= entry_count
                ),
                aa_counts_json TEXT NOT NULL,
                target_id_fingerprint TEXT,
                target_description_fingerprint TEXT,
                target_content_fingerprint TEXT,
                upper_cased_entries INTEGER NOT NULL CHECK (upper_cased_entries >= 0),
                stop_stripped_entries INTEGER NOT NULL CHECK (stop_stripped_entries >= 0),
                illegal_residue_entries INTEGER NOT NULL CHECK (illegal_residue_entries >= 0),
                illegal_residues_json TEXT NOT NULL,
                empty_sequence_entries INTEGER NOT NULL CHECK (empty_sequence_entries >= 0),
                bare_identifier_entries INTEGER NOT NULL CHECK (bare_identifier_entries >= 0),
                id_namespaces_json TEXT NOT NULL,
                CHECK (
                    (detail_level = 'full'
                        AND distinct_target_ids IS NOT NULL
                        AND distinct_target_sequences IS NOT NULL
                        AND distinct_target_descriptions IS NOT NULL
                        AND duplicate_target_id_occurrences IS NOT NULL
                        AND conflicting_target_ids IS NOT NULL
                        AND repeated_target_sequences IS NOT NULL
                        AND target_id_fingerprint IS NOT NULL
                        AND target_description_fingerprint IS NOT NULL
                        AND target_content_fingerprint IS NOT NULL
                        AND aa_sample_size = entry_count)
                    OR
                    (detail_level = 'metadata_only'
                        AND distinct_target_ids IS NULL
                        AND distinct_target_sequences IS NULL
                        AND distinct_target_descriptions IS NULL
                        AND duplicate_target_id_occurrences IS NULL
                        AND conflicting_target_ids IS NULL
                        AND repeated_target_sequences IS NULL
                        AND target_id_fingerprint IS NULL
                        AND target_description_fingerprint IS NULL
                        AND target_content_fingerprint IS NULL)
                )
            )
        """,
    "database_kind_stats": """
        CREATE TABLE IF NOT EXISTS database_kind_stats (
                database_id INTEGER NOT NULL REFERENCES databases(id) ON DELETE CASCADE,
                kind TEXT NOT NULL CHECK (kind IN ('target', 'decoy', 'contaminant', 'entrapment', 'sentinel')),
                entry_count INTEGER NOT NULL CHECK (entry_count >= 0),
                distinct_ids INTEGER CHECK (distinct_ids IS NULL OR distinct_ids >= 0),
                distinct_sequences INTEGER CHECK (distinct_sequences IS NULL OR distinct_sequences >= 0),
                distinct_descriptions INTEGER CHECK (distinct_descriptions IS NULL OR distinct_descriptions >= 0),
                distinct_pairs INTEGER CHECK (distinct_pairs IS NULL OR distinct_pairs >= 0),
                duplicate_id_occurrences INTEGER
                    CHECK (duplicate_id_occurrences IS NULL OR duplicate_id_occurrences >= 0),
                conflicting_ids INTEGER CHECK (conflicting_ids IS NULL OR conflicting_ids >= 0),
                repeated_sequences INTEGER CHECK (repeated_sequences IS NULL OR repeated_sequences >= 0),
                length_min INTEGER NOT NULL CHECK (length_min >= 0),
                length_q1 REAL NOT NULL CHECK (length_q1 >= 0),
                length_median REAL NOT NULL CHECK (length_median >= 0),
                length_mean REAL NOT NULL CHECK (length_mean >= 0),
                length_q3 REAL NOT NULL CHECK (length_q3 >= 0),
                length_max INTEGER NOT NULL CHECK (length_max >= 0),
                total_residues INTEGER NOT NULL CHECK (total_residues >= 0),
                aa_sample_size INTEGER NOT NULL CHECK (
                    aa_sample_size >= 0 AND aa_sample_size <= entry_count
                ),
                aa_counts_json TEXT NOT NULL,
                id_fingerprint TEXT,
                description_fingerprint TEXT,
                content_fingerprint TEXT,
                PRIMARY KEY (database_id, kind)
            ) WITHOUT ROWID
        """,
    "entries": """
        CREATE TABLE IF NOT EXISTS entries (
                database_id INTEGER NOT NULL REFERENCES databases(id) ON DELETE CASCADE,
                ordinal INTEGER NOT NULL,
                sequence_id TEXT NOT NULL COLLATE BINARY,
                kind TEXT NOT NULL CHECK (kind IN ('target', 'decoy', 'contaminant', 'entrapment', 'sentinel')),
                contaminant_group TEXT,
                sequence_length INTEGER NOT NULL CHECK (sequence_length >= 0),
                sequence_hash BLOB NOT NULL CHECK (length(sequence_hash) = 16),
                description_hash BLOB CHECK (description_hash IS NULL OR length(description_hash) = 16),
                PRIMARY KEY (database_id, ordinal),
                CHECK (
                    (kind = 'contaminant' AND contaminant_group IS NOT NULL)
                    OR (kind != 'contaminant' AND contaminant_group IS NULL)
                )
            ) WITHOUT ROWID
        """,
    "database_pair_stats": """
        CREATE TABLE IF NOT EXISTS database_pair_stats (
                database_id_low INTEGER NOT NULL
                    REFERENCES databases(id) ON DELETE CASCADE,
                database_id_high INTEGER NOT NULL
                    REFERENCES databases(id) ON DELETE CASCADE,
                kind TEXT NOT NULL CHECK (kind IN ('target', 'contaminant')),
                shared_ids INTEGER NOT NULL CHECK (shared_ids >= 0),
                shared_sequence_checksums INTEGER NOT NULL
                    CHECK (shared_sequence_checksums >= 0),
                shared_descriptions INTEGER NOT NULL CHECK (shared_descriptions >= 0),
                shared_exact_pairs INTEGER NOT NULL CHECK (shared_exact_pairs >= 0),
                matching_shared_ids INTEGER NOT NULL CHECK (matching_shared_ids >= 0),
                PRIMARY KEY (database_id_low, database_id_high, kind),
                CHECK (database_id_low < database_id_high),
                CHECK (matching_shared_ids <= shared_ids),
                CHECK (matching_shared_ids <= shared_exact_pairs)
            ) WITHOUT ROWID
        """,
}
"""DDL per registry table, keyed by name so a partial copy can select tables."""


KIND_STATS_FINGERPRINT_INDEX = "CREATE INDEX IF NOT EXISTS database_kind_stats_kind_id_fingerprint ON database_kind_stats(kind, id_fingerprint)"
PAIR_STATS_HIGH_INDEX = (
    "CREATE INDEX IF NOT EXISTS database_pair_stats_high "
    "ON database_pair_stats(database_id_high, database_id_low, kind)"
)

ENTRY_LOOKUP_INDEXES: tuple[str, ...] = (
    "CREATE INDEX IF NOT EXISTS entries_kind_db_id ON entries(kind, database_id, sequence_id)",
    "CREATE INDEX IF NOT EXISTS entries_kind_id_db ON entries(kind, sequence_id, database_id)",
    "CREATE INDEX IF NOT EXISTS entries_kind_db_sequence ON entries(kind, database_id, sequence_hash)",
    "CREATE INDEX IF NOT EXISTS entries_kind_sequence_db ON entries(kind, sequence_hash, database_id)",
    "CREATE INDEX IF NOT EXISTS entries_kind_db_description ON entries(kind, database_id, description_hash) "
    "WHERE description_hash IS NOT NULL",
    "CREATE INDEX IF NOT EXISTS entries_kind_description_db ON entries(kind, description_hash, database_id) "
    "WHERE description_hash IS NOT NULL",
    KIND_STATS_FINGERPRINT_INDEX,
)

PAIR_LOOKUP_INDEXES: tuple[str, ...] = (PAIR_STATS_HIGH_INDEX,)

STATS_ONLY_INDEXES: tuple[str, ...] = (KIND_STATS_FINGERPRINT_INDEX, PAIR_STATS_HIGH_INDEX)
"""The indexes belonging to the tables a stats-only copy carries."""


CHILD_TABLES: tuple[tuple[str, str], ...] = (
    ("entries", "database_id"),
    ("database_kind_stats", "database_id"),
    ("database_pair_stats", "database_id_low"),
    ("database_pair_stats", "database_id_high"),
)

UNIQUE_KEYS: tuple[tuple[str, str], ...] = (
    ("entries", "database_id, ordinal"),
    ("database_kind_stats", "database_id, kind"),
    ("database_pair_stats", "database_id_low, database_id_high, kind"),
    ("registry_meta", "key"),
)
