class Conservation
  include DataMapper::Resource
  
  property :id, Serial, :key=>true
  property :alignment_id, Integer, :required => true
  property :alignment_name, String, :required => true, :key=>true
  property :conservation_results, Text, :required => true
  property :deleted_at, ParanoidDateTime
  
  def results_array
    self.conservation_results.split()[1..-1]
  end
  
  def type
    self.conservation_results.split()[0]
  end
  
end
