def get_library(cursor, sample, type):
    """Get the library ID for a sample"""
    query = "SELECT id FROM sample WHERE sample_id = %s"
    cursor.execute(query, sample)
    sample_id = cursor.fetchone()[0]
    
    query = "SELECT id FROM library WHERE sample_id = %s AND library_type = %s"
    cursor.execute(query, (sample_id, type))
    return cursor.fetchone()[0]

def addEvent(cursor, library_id, event_type):
    """Add a generic event to the database"""
    query = "INSERT INTO event (library_id, type) VALUES (%s, %s)"
    cursor.execute(query, (library_id, event_type))
    return cursor.lastrowid

LIBRARY_TYPES = ["exome", "genome", "RNA-seq", "enriched", "meta"]
MUTATION_TYPES = ["germline", "somatic", "artifact", "unknown"]
