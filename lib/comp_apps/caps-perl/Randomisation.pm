package Randomisation;
require Exporter;
@ISA = qw(Exporter);
	
my(@codigo_genetico) = ('F','L','I','M','V','S','P','T','A','Y','H','Q','N','K','D','E','C','W','R','G');
sub mutate
{
	my($seq, $mutations) = @_;
	my($position) = random_position($seq);
	my($aa) = random_amino(@codigo_genetico);
	substr($seq, $position, 1, $aa);
	return $seq;
}
sub random_position
{
	my($sequence) = @_;
	return (int(rand(length($sequence) - 1)));
}
sub random_amino
{
	my($AA) = $codigo_genetico[(int(rand(scalar(@codigo_genetico) - 1)))];
	return $AA;
}